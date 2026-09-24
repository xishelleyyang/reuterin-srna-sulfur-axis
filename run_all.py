#!/usr/bin/env python
"""Orchestrate the full reuterin sRNA/sulfur-axis reanalysis, end to end.

    python run_all.py --data-dir data --out-dir results --stage all
    python run_all.py --data-dir data --out-dir results --stage transcriptomes --pilot
    python run_all.py --data-dir data --out-dir results --stage figures

Stages (run in this order for --stage all):
  transcriptomes   Rscript src/recompute_statistics.R  (microarray + RNA-seq DE, gene-set roast tests)
  conservation     src/conservation_analysis.py + src/typhimurium_duplex_mfe.py
  tables           src/build_tables.py                 (assembles Tables S1-S7 from the two stages above)
  figures          src/make_figures.py + src/make_fig5_conservation.py
  publish          copy the 7 supplementary tables and 5 figures to the filenames actually
                    cited in the manuscript (source scripts keep their own internal "_v4"
                    tag; this stage does the one mechanical rename, transparently, so nothing
                    inside a verified script needs to be edited to match a later manuscript
                    version bump)
  validate         re-check output files exist, are non-empty, and spot-check a few
                    manuscript-reported numbers (see tests/test_pipeline_smoke.py for the
                    full set)

--pilot routes every stage's output under results/pilot/ instead of results/formal/, and
runs the R gene-set rotation tests at NROT=199 (vs. 9999) for a fast sanity pass. Pilot
output is never mixed into the formal tree and the validate/publish stages both refuse to
treat a pilot run as submission-ready.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "src"

# (internal "_v4" name written by build_tables.py / make_figures.py, manuscript-cited name)
PUBLISHED_TABLES = [
    ("Table_S1_CLASH_interactions_v4.csv", "Table_S1_CLASH_interactions_v5.csv"),
    ("Table_S2_gene_set_statistics_v4.csv", "Table_S2_gene_set_statistics_v5.csv"),
    ("Table_S3_conservation_v4.csv", "Table_S3_conservation_v5.csv"),
    ("Table_S4_Reuterin_DE_v4.csv", "Table_S4_Reuterin_DE_v5.csv"),
    ("Table_S5_network_sensitivity_v4.csv", "Table_S5_network_sensitivity_v5.csv"),
    ("Table_S6_gene_resampling_v4.csv", "Table_S6_gene_resampling_v5.csv"),
    ("Table_S7_block_sensitivity_v4.csv", "Table_S7_block_sensitivity_v5.csv"),
]
PUBLISHED_FIGURE_STEMS = [
    ("Fig1_reuterin_program_v4", "Fig1_reuterin_program_v5"),
    ("Fig2_isc_operon_v4", "Fig2_isc_operon_v5"),
    ("Fig3_candidate_network_v4", "Fig3_candidate_network_v5"),
    ("Fig4_crossstress_v4", "Fig4_crossstress_v5"),
    ("Fig5_conservation_v4", "Fig5_conservation_v5"),
]


def stage_dirs(out_dir: Path, pilot: bool) -> dict[str, Path]:
    root = out_dir / ("pilot" if pilot else "formal")
    final = root / "final"
    return {
        "root": root,
        "stats": root / "statistics",
        "compgen": root / "compgen",
        "final": final,
        "supplementary": final / "supplementary",
        "figures": final / "figures",
        "published": root / "published",
    }


def run(cmd: list[str]) -> None:
    print("+", " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True)


def stage_transcriptomes(data_dir: Path, dirs: dict, pilot: bool) -> None:
    cmd = ["Rscript", str(SRC / "recompute_statistics.R"), str(data_dir), str(dirs["stats"])]
    if pilot:
        cmd.append("--pilot")
    run(cmd)
    if pilot:
        # recompute_statistics.R --pilot uses a single seed (SEEDS <- 42L), so it never
        # writes genesets_seed_check.csv (that file only exists when a second seed ran).
        # build_tables.py unconditionally reads it to build the seed-stability column.
        # This asymmetry never affects any manuscript number -- every reported q-value
        # comes from the formal (non-pilot, dual-seed) run, which always writes the file
        # normally -- but it does mean pilot output cannot be handed to build_tables.py
        # as-is. Rather than edit either verified script, stand in a placeholder here so
        # pilot mode can smoke-test the full chain; the resulting "second_seed_q_up" /
        # "seed_threshold_stable" columns in pilot output are trivially self-identical
        # and carry no scientific meaning -- only the formal run's are real.
        seed_check = dirs["stats"] / "genesets_seed_check.csv"
        if not seed_check.exists():
            shutil.copy(dirs["stats"] / "genesets_primary.csv", seed_check)
            print(f"[pilot] wrote placeholder {seed_check.name} "
                  "(single-seed pilot run has no independent second seed; "
                  "not meaningful, smoke-test only)")


def stage_conservation(data_dir: Path, dirs: dict) -> None:
    run([sys.executable, str(SRC / "conservation_analysis.py"),
         "--data", str(data_dir), "--out", str(dirs["compgen"])])
    run([sys.executable, str(SRC / "typhimurium_duplex_mfe.py"),
         "--data", str(data_dir), "--compgen-out", str(dirs["compgen"]), "--out", str(dirs["compgen"])])


def stage_tables(data_dir: Path, dirs: dict) -> None:
    dirs["final"].mkdir(parents=True, exist_ok=True)
    run([sys.executable, str(SRC / "build_tables.py"),
         "--data", str(data_dir), "--stats", str(dirs["stats"]),
         "--compgen", str(dirs["compgen"]), "--out", str(dirs["final"])])


def stage_figures(data_dir: Path, dirs: dict) -> None:
    run([sys.executable, str(SRC / "make_figures.py"),
         "--stats", str(dirs["stats"]), "--final", str(dirs["final"]),
         "--data", str(data_dir), "--fig-out", str(dirs["figures"])])
    run([sys.executable, str(SRC / "make_fig5_conservation.py"),
         "--compgen-out", str(dirs["compgen"]), "--fig-out", str(dirs["figures"])])


def stage_publish(dirs: dict, pilot: bool) -> None:
    dirs["published"].mkdir(parents=True, exist_ok=True)
    n = 0
    for src_name, dst_name in PUBLISHED_TABLES:
        src = dirs["supplementary"] / src_name
        if not src.exists():
            raise FileNotFoundError(f"expected table not found: {src} (run stage=tables first)")
        shutil.copy(src, dirs["published"] / dst_name)
        n += 1
    for src_stem, dst_stem in PUBLISHED_FIGURE_STEMS:
        for ext in ("png", "svg"):
            src = dirs["figures"] / f"{src_stem}.{ext}"
            if not src.exists():
                raise FileNotFoundError(f"expected figure not found: {src} (run stage=figures first)")
            shutil.copy(src, dirs["published"] / f"{dst_stem}.{ext}")
            n += 1
    extra = dirs["compgen"] / "Table_S3b_duplex_MFE.csv"
    if extra.exists():
        shutil.copy(extra, dirs["published"] / extra.name)
        n += 1
    print(f"published {n} file(s) to {dirs['published']}")
    if pilot:
        (dirs["published"] / "PILOT_DO_NOT_SUBMIT.txt").write_text(
            "This directory was produced with --pilot (NROT=199 rotation tests). "
            "It is a fast sanity check, not the submission-grade result. Re-run without "
            "--pilot before using any number or figure from here.\n")
        print("PILOT run: wrote PILOT_DO_NOT_SUBMIT.txt -- do not submit these files.")


def stage_validate(dirs: dict, pilot: bool) -> None:
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "test_pipeline_smoke", HERE / "tests" / "test_pipeline_smoke.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.run_all_checks(dirs["root"], pilot=pilot)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", default="data")
    p.add_argument("--out-dir", default="results")
    p.add_argument("--stage", default="all",
                   choices=["transcriptomes", "conservation", "tables", "figures", "publish", "validate", "all"])
    p.add_argument("--pilot", action="store_true",
                   help="fast sanity pass (NROT=199); routed to results/pilot/, never results/formal/")
    a = p.parse_args()

    data_dir = Path(a.data_dir).resolve()
    out_dir = Path(a.out_dir).resolve()
    dirs = stage_dirs(out_dir, a.pilot)

    order = ["transcriptomes", "conservation", "tables", "figures", "publish", "validate"]
    stages = order if a.stage == "all" else [a.stage]

    for stage in stages:
        print(f"\n=== stage: {stage} ({'pilot' if a.pilot else 'formal'}) ===")
        if stage == "transcriptomes":
            stage_transcriptomes(data_dir, dirs, a.pilot)
        elif stage == "conservation":
            stage_conservation(data_dir, dirs)
        elif stage == "tables":
            stage_tables(data_dir, dirs)
        elif stage == "figures":
            stage_figures(data_dir, dirs)
        elif stage == "publish":
            stage_publish(dirs, a.pilot)
        elif stage == "validate":
            stage_validate(dirs, a.pilot)

    print(f"\nDone. Output root: {dirs['root']}")


if __name__ == "__main__":
    main()
