"""Regression guard for the curated gene-set definitions.

This does not re-derive biology; it protects against the exact class of bug
found during the manuscript audit: geneset_unified_v3.R (now retired, not
part of this repository) used the invalid symbol "gorA" instead of "gor",
and did not exclude the Hpx-strain-deleted genes (katE, katG, ahpC) from
gene-set scoring. Both defects silently changed a reported log2FC/q-value.
These checks parse the gene sets directly out of the shipped source files
(not a second hand-typed copy) and fail loudly if either defect reappears,
or if any set's membership changes without a matching manuscript update.

Run directly: python tests/test_gene_sets.py
Or via pytest: pytest tests/test_gene_sets.py
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src"
sys.path.insert(0, str(SRC))

# Fingerprint of the 8 curated gene sets as they must read in recompute_statistics.R,
# i.e. the set actually used for every roast()/gene-set number quoted in the manuscript.
EXPECTED_SETS_SHA256 = "6ddd9517ea7a92f8df24c05bd24119fc0e95c1366b8f6109a5f760f279d72d9e"
HPX_DELETED_EXPECTED = {"katE", "katG", "ahpC"}


def parse_r_sets(r_source: str) -> dict[str, list[str]]:
    block = re.search(r"SETS\s*<-\s*list\((.*?)\)\n(?=HPX_DELETED)", r_source, re.S)
    assert block, "could not locate the SETS <- list(...) block in recompute_statistics.R"
    body = block.group(1)
    sets = {}
    for name, members in re.findall(r'"([^"]+)"\s*=\s*c\(([^)]*)\)', body):
        sets[name] = re.findall(r'"([^"]+)"', members)
    return sets


def parse_hpx_deleted(r_source: str) -> set[str]:
    m = re.search(r'HPX_DELETED\s*<-\s*c\(([^)]*)\)', r_source)
    assert m, "could not locate HPX_DELETED <- c(...) in recompute_statistics.R"
    return set(re.findall(r'"([^"]+)"', m.group(1)))


def test_sets_fingerprint_unchanged():
    r_source = (SRC / "recompute_statistics.R").read_text()
    sets = parse_r_sets(r_source)
    assert len(sets) == 8, f"expected 8 curated gene sets, found {len(sets)}: {sorted(sets)}"
    fp = hashlib.sha256(json.dumps(sets, sort_keys=True).encode()).hexdigest()
    assert fp == EXPECTED_SETS_SHA256, (
        f"gene-set membership changed (fingerprint {fp} != expected {EXPECTED_SETS_SHA256}). "
        "If this is an intentional, manuscript-approved change, update EXPECTED_SETS_SHA256 "
        "deliberately -- do not silently relax this check.")


def test_no_stale_gorA_alias():
    r_source = (SRC / "recompute_statistics.R").read_text()
    sets = parse_r_sets(r_source)
    all_symbols = {s for members in sets.values() for s in members}
    assert "gorA" not in all_symbols, "stale invalid symbol 'gorA' reappeared (canonical symbol is 'gor', b3500)"
    assert "gor" in sets["OxyR regulon"] and "gor" in sets["Glutathione/glyoxalase"]


def test_hpx_deleted_genes_excluded_from_scoring():
    r_source = (SRC / "recompute_statistics.R").read_text()
    hpx_deleted = parse_hpx_deleted(r_source)
    assert hpx_deleted == HPX_DELETED_EXPECTED
    # Confirm the exclusion logic is actually wired to the Hpx contrast, not just declared.
    assert re.search(r'exc\s*<-\s*if\(ct=="Hpx"\)\s*HPX_DELETED', r_source), \
        "HPX_DELETED no longer appears to be excluded specifically for the Hpx contrast"
    sets = parse_r_sets(r_source)
    affected = {s for s in sets["OxyR regulon"] if s in hpx_deleted}
    assert affected == {"ahpC", "katG"}, (
        "expected exactly {ahpC, katG} of the 3 Hpx-deleted genes to sit inside the OxyR "
        f"regulon set as curated; found {affected}")


def test_pathways_dict_matches_build_tables():
    import build_tables  # noqa: E402  (path inserted above)
    pathways = build_tables.PATHWAYS
    sulfur = set(sum(pathways.values(), []))
    fes = set(pathways["ISC"] + pathways["SUF"] + pathways["FeS/iron"])
    assert len(sulfur) == 86, f"sulfur universe size drifted: {len(sulfur)} (expected 86)"
    assert len(fes) == 29, f"Fe-S/iron universe size drifted: {len(fes)} (expected 29)"
    assert "gor" in pathways["Glutathione"] and "gorA" not in sulfur


def run_all() -> None:
    tests = [v for k, v in globals().items() if k.startswith("test_") and callable(v)]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS  {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {t.__name__}: {e}")
    if failed:
        raise SystemExit(f"{failed}/{len(tests)} gene-set regression checks failed")
    print(f"All {len(tests)} gene-set regression checks passed.")


if __name__ == "__main__":
    run_all()
