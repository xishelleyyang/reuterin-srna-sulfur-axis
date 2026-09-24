#!/usr/bin/env python
"""RyhB-iscS duplex free energy for MG1655 and S. Typhimurium LT2 (Discussion).

The manuscript Discussion states: "the predicted RyhB-iscS duplex for
S. Typhimurium LT2 remains highly stable (minimum free energy -76.0 versus
-54.7 kcal/mol for MG1655)". conservation_analysis.py computes the MG1655
side of that comparison directly (it is the BLAST query), but does not by
itself emit a per-target-genome duplex value -- no script in the project's
history computed the LT2 side as a saved, reusable step; it existed only as
a one-off validation calculation. This script closes that gap using the
*same* ViennaRNA RNA.duplexfold call already used for MG1655, applied to
S. Typhimurium LT2's own genomic ryhB-1 and iscS 5'UTR sequences. Those
sequences are recovered by stripping alignment gap characters ('-') from the
LT2 entry in alignments_v3.json (conservation_analysis.py's output) -- i.e.
by gap-stripping the alignment BLAST/PairwiseAligner already produced, not
by re-deriving a new alignment or applying any new method.

Scope: this script intentionally reproduces only the two values quoted in
the manuscript text (MG1655 reference and S. Typhimurium LT2). alignments_v3.json
contains the same gap-containing subject alignment for all 13 target genomes,
so the identical duplexfold-on-gap-stripped-sequence procedure could be
extended to any of the other 11 genomes; that extension is deliberately not
performed here because it is not a value reported anywhere in the manuscript,
and shipping unreported numbers under a manuscript-facing table risks being
mistaken for a published claim. Use --target-genome to compute a different
single genome for exploratory purposes if needed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import RNA

from mg1655_queries import load_mg1655_queries

DEFAULT_TARGET = "S.Typhimurium_LT2"


def duplex_mfe_for_target(compgen_out: Path, genome: str) -> dict:
    alignments = json.loads((compgen_out / "alignments_v3.json").read_text())
    if genome not in alignments:
        raise KeyError(
            f"{genome!r} not found in alignments_v3.json (available: "
            f"{sorted(alignments)}). Run conservation_analysis.py first.")
    a = alignments[genome]
    ryhB_g = a["ryhB_s"].replace("-", "")
    iscS_g = a["iscS_s"].replace("-", "")
    dup = RNA.duplexfold(ryhB_g, iscS_g)
    return dict(genome=genome, ryhB_gapstripped_nt=len(ryhB_g), iscS_gapstripped_nt=len(iscS_g),
                mfe_kcal_mol=round(dup.energy, 1), structure=dup.structure)


def main(data: str | Path, compgen_out: str | Path, out: str | Path, target_genome: str) -> pd.DataFrame:
    data, compgen_out, out = Path(data), Path(compgen_out), Path(out)
    out.mkdir(parents=True, exist_ok=True)

    ryhB, iscS_utr = load_mg1655_queries(data)
    ref = RNA.duplexfold(ryhB, iscS_utr)
    rows = [dict(genome="MG1655 (reference query)", ryhB_gapstripped_nt=len(ryhB),
                 iscS_gapstripped_nt=len(iscS_utr), mfe_kcal_mol=round(ref.energy, 1),
                 structure=ref.structure)]
    rows.append(duplex_mfe_for_target(compgen_out, target_genome))

    df = pd.DataFrame(rows)
    df.to_csv(out / "Table_S3b_duplex_MFE.csv", index=False)
    print(df[["genome", "ryhB_gapstripped_nt", "iscS_gapstripped_nt", "mfe_kcal_mol"]].to_string(index=False))
    return df


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", required=True, help="directory with mg1655.fna.gz")
    p.add_argument("--compgen-out", required=True, help="output directory from conservation_analysis.py")
    p.add_argument("--out", required=True, help="directory to write Table_S3b_duplex_MFE.csv")
    p.add_argument("--target-genome", default=DEFAULT_TARGET,
                   help=f"genome key in alignments_v3.json to compare against MG1655 (default: {DEFAULT_TARGET})")
    a = p.parse_args()
    main(a.data, a.compgen_out, a.out, a.target_genome)
