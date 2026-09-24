"""Shared MG1655 query-sequence extraction for the conservation module.

Both conservation_analysis.py (all 13 genomes) and typhimurium_duplex_mfe.py
(duplex free-energy for the two genomes discussed in the manuscript text) need
the identical MG1655 ryhB and iscS 5'UTR+TIR query sequences. Factored out here
so the two scripts cannot silently drift onto different coordinates.

Coordinates are unchanged from the original module (0-based, half-open, on the
MG1655 chromosome, minus strand for both features):
  ryhB:     [3580921:3581016)  -> 95 nt
  iscS_utr: [2661471:2661681)  -> 210 nt, position 150 (0-based) is the iscS ATG
"""
from __future__ import annotations

import gzip
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

RYHB_SLICE = (3580921, 3581016)
ISCS_UTR_SLICE = (2661471, 2661681)
ISCS_ATG_OFFSET = 150  # 0-based position of the iscS start codon within iscS_utr


def load_mg1655_genome(data_dir: str | Path, filename: str = "mg1655.fna.gz") -> str:
    path = Path(data_dir) / filename
    rec = next(SeqIO.parse(gzip.open(path, "rt"), "fasta"))
    return str(rec.seq).upper()


def load_mg1655_queries(data_dir: str | Path, filename: str = "mg1655.fna.gz") -> tuple[str, str]:
    """Return (ryhB, iscS_utr) exactly as used by the frozen analysis."""
    genome = load_mg1655_genome(data_dir, filename)
    lo, hi = RYHB_SLICE
    ryhB = str(Seq(genome[lo:hi]).reverse_complement())
    lo, hi = ISCS_UTR_SLICE
    iscS_utr = str(Seq(genome[lo:hi]).reverse_complement())
    assert iscS_utr[ISCS_ATG_OFFSET:ISCS_ATG_OFFSET + 3] == "ATG", "start codon check failed"
    assert len(ryhB) == 95 and len(iscS_utr) == 210
    return ryhB, iscS_utr
