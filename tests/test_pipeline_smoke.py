"""End-to-end smoke test: did a run_all.py execution actually produce the expected
files, and do a handful of manuscript-reported numbers still come out right?

This is intentionally NOT a full re-derivation (that is run_all.py itself, and the
audit report's independent reproduction). It is the fast, automatable check that a
--pilot or formal run really executed the whole chain rather than silently stopping
partway, and that a handful of headline numbers have not drifted.

Usage:
    pytest tests/test_pipeline_smoke.py --results-root results/formal
    python tests/test_pipeline_smoke.py results/formal --pilot     # standalone, or called by run_all.py

Numeric tolerance follows the plan's regression convention: rtol=1e-7 for
full-precision values recomputed here; manuscript-rounded values (e.g. "-54.7
kcal/mol") are checked at the precision they are actually reported in text.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd

REQUIRED_FILES_FORMAL_ONLY = [
    # roast() at NROT=9999 is what the manuscript's q-values were computed with;
    # a --pilot run (NROT=199) intentionally will not match these to full precision.
]


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise AssertionError(f"missing expected output file: {path}")
    if path.stat().st_size == 0:
        raise AssertionError(f"output file is empty: {path}")
    return pd.read_csv(path)


def check_files_exist(root: Path) -> list[str]:
    problems = []
    expected = [
        root / "statistics" / "DE_Reuterin.csv",
        root / "statistics" / "DE_HOCl.csv",
        root / "statistics" / "DE_Ferrate.csv",
        root / "statistics" / "DE_Hpx.csv",
        root / "statistics" / "genesets_primary.csv",
        root / "compgen" / "Table_S3_conservation_v3.csv",
        root / "compgen" / "alignments_v3.json",
        root / "compgen" / "Table_S3b_duplex_MFE.csv",
        root / "final" / "supplementary" / "Table_S1_CLASH_interactions_v4.csv",
        root / "final" / "supplementary" / "Table_S4_Reuterin_DE_v4.csv",
        root / "final" / "audit" / "result_summary_v4.json",
        root / "final" / "figures" / "Fig1_reuterin_program_v4.png",
        root / "final" / "figures" / "Fig5_conservation_v4.png",
        root / "published" / "Table_S1_CLASH_interactions_v5.csv",
        root / "published" / "Fig5_conservation_v5.png",
    ]
    for f in expected:
        if not f.exists():
            problems.append(f"missing: {f}")
        elif f.stat().st_size == 0:
            problems.append(f"empty: {f}")
    return problems


def check_de_reuterin_counts(root: Path) -> list[str]:
    """4,462 genes; 71 at 5% FDR; 117 at 10% FDR (Results paragraph 1)."""
    problems = []
    de = _read_csv(root / "statistics" / "DE_Reuterin.csv")
    n = len(de)
    fdr05 = int((de["adj.P.Val"] < 0.05).sum())
    fdr10 = int((de["adj.P.Val"] < 0.10).sum())
    if n != 4462:
        problems.append(f"DE_Reuterin.csv row count = {n}, expected 4462")
    if fdr05 != 71:
        problems.append(f"DE_Reuterin FDR<0.05 count = {fdr05}, expected 71")
    if fdr10 != 117:
        problems.append(f"DE_Reuterin FDR<0.10 count = {fdr10}, expected 117")
    return problems


def check_geneset_roast(root: Path, pilot: bool) -> list[str]:
    """Electrophile detox +3.07 q=0.008 and ISC +1.59 q=0.012 in the Reuterin
    contrast (Results paragraph 1). q-value match requires the formal NROT=9999
    run; a pilot run at NROT=199 is checked for direction and rough magnitude only."""
    problems = []
    p = _read_csv(root / "statistics" / "genesets_primary.csv")
    reu = p[p.contrast == "Reuterin"].set_index("gene_set")
    checks = [("Electrophile detoxification", 3.07, 0.008), ("ISC operon", 1.59, 0.012)]
    for gene_set, exp_fc, exp_q in checks:
        if gene_set not in reu.index:
            problems.append(f"gene set {gene_set!r} missing from genesets_primary.csv")
            continue
        row = reu.loc[gene_set]
        if not math.isclose(row.raw_mean_logFC, exp_fc, abs_tol=0.02):
            problems.append(f"{gene_set} raw_mean_logFC={row.raw_mean_logFC:.3f}, expected ~{exp_fc}")
        if not pilot and not math.isclose(row.roast_q_up, exp_q, abs_tol=0.01):
            problems.append(f"{gene_set} roast_q_up={row.roast_q_up:.4f}, expected ~{exp_q} (formal NROT run)")
    return problems


def check_clash_enrichment(root: Path) -> list[str]:
    """Fe-S/iron subsystem OR=4.44, P=8.34e-5 (Results paragraph on the CLASH network)."""
    problems = []
    t5 = _read_csv(root / "final" / "supplementary" / "Table_S5_network_sensitivity_v4.csv")
    row = t5[(t5.analysis == "Primary") & (t5.subsystem == "FeS/iron")]
    if len(row) != 1:
        problems.append(f"expected exactly 1 Primary/FeS-iron row in Table S5, found {len(row)}")
        return problems
    row = row.iloc[0]
    if not math.isclose(row.OR, 4.44, abs_tol=0.02):
        problems.append(f"Fe-S/iron OR={row.OR:.3f}, expected ~4.44")
    if not math.isclose(row.P, 8.344645e-05, rel_tol=1e-3):
        problems.append(f"Fe-S/iron P={row.P:.3e}, expected ~8.3e-5")
    return problems


def check_duplex_mfe(root: Path) -> list[str]:
    """MG1655 reference -54.7 kcal/mol; S. Typhimurium LT2 -76.0 kcal/mol (Discussion)."""
    problems = []
    df = _read_csv(root / "compgen" / "Table_S3b_duplex_MFE.csv").set_index("genome")
    expected = {"MG1655 (reference query)": -54.7, "S.Typhimurium_LT2": -76.0}
    for genome, exp_mfe in expected.items():
        if genome not in df.index:
            problems.append(f"{genome!r} missing from Table_S3b_duplex_MFE.csv")
            continue
        got = df.loc[genome, "mfe_kcal_mol"]
        if not math.isclose(got, exp_mfe, abs_tol=0.05):
            problems.append(f"{genome} MFE={got}, expected {exp_mfe} kcal/mol")
    return problems


def run_all_checks(root: Path, pilot: bool = False) -> None:
    root = Path(root)
    problems: list[str] = []
    problems += check_files_exist(root)
    if not problems:  # only chase numbers once the files we need actually exist
        problems += check_de_reuterin_counts(root)
        problems += check_geneset_roast(root, pilot)
        problems += check_clash_enrichment(root)
        problems += check_duplex_mfe(root)

    if problems:
        msg = "\n".join(f"  - {p}" for p in problems)
        raise AssertionError(f"{len(problems)} smoke-test check(s) failed under {root}:\n{msg}")
    tag = "PILOT" if pilot else "FORMAL"
    print(f"[{tag}] all smoke-test checks passed under {root}")


if __name__ == "__main__":
    argv = [a for a in sys.argv[1:] if a != "--pilot"]
    pilot = "--pilot" in sys.argv
    if not argv:
        raise SystemExit("usage: python test_pipeline_smoke.py RESULTS_ROOT [--pilot]")
    run_all_checks(Path(argv[0]), pilot=pilot)
