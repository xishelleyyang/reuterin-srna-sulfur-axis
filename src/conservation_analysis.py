#!/usr/bin/env python
"""Module 4: conservation of ryhB and the iscS 5'UTR across 13 target genomes,
relative to the MG1655 reference (Fig. 5; Table S3).

Method notes (why this is not a naive single best-hit BLAST report):
 (1) paralog-aware ryhB hit reporting -- every Salmonella genome carries a
     second, partial ryhB copy in addition to the syntenic full-length locus,
     and both are reported rather than only the top hit;
 (2) explicit gap handling -- identity is reported three ways (per query
     length / per aligned column / raw gap count) so that a deletion is never
     silently counted as either a match or a mismatch;
 (3) a symmetric 15-nt sliding window with a >=70% coverage requirement is
     used for the smoothed per-position conservation trace, avoiding the
     edge-decay artifact of a causal/one-sided window;
 (4) pairing-window identity is reported per species group rather than as a
     single pooled percentage.

Inputs (under --data):
  mg1655.fna.gz                  MG1655 RefSeq genome (GCF_000005845.2)
  genomes/<name>.fna.gz           one gzipped FASTA per target genome

Outputs (under --out):
  queries.fasta                   the two MG1655 query sequences, for the record
  hits_all_<genome>.txt            raw blastn hits per genome
  Table_S3_conservation_v3.csv     per-genome summary (consumed by build_tables.py)
  alignments_v3.json               per-genome aligned query/subject strings
                                    (consumed by build_tables.py and
                                    typhimurium_duplex_mfe.py)
  conservation_paralogs.json        per-genome paralog-resolved hit detail
  iscS_perposition_conservation_v3.csv  position-resolved identity trace

The output filenames intentionally keep their internal "_v3" tag: they are
consumed programmatically by build_tables.py and typhimurium_duplex_mfe.py,
not cited anywhere in the manuscript text, so there is no reason to rename
them and every reason to keep the interface between scripts unchanged from
the verified pipeline.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import RNA
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

from mg1655_queries import load_mg1655_queries

COLS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore"]


class NpEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o) if not np.isnan(o) else None
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)


def run_blast(data: Path, out: Path, genomes: list[str]) -> None:
    outfmt = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
    for g in genomes:
        fna = out / f".{g}.fna"
        subprocess.run(f"gunzip -c {data}/genomes/{g}.fna.gz > {fna}", shell=True, check=True)
        subprocess.run(
            f"makeblastdb -in {fna} -dbtype nucl -out {out}/.{g} -logfile /dev/null",
            shell=True, check=True)
        subprocess.run(
            f'blastn -query {data}/queries.fasta -db {out}/.{g} -outfmt "{outfmt}" '
            f'-word_size 7 -evalue 1e-3 -dust no > {out}/hits_all_{g}.txt',
            shell=True, check=True)
        subprocess.run(f"rm -f {fna} {out}/.{g}.n*", shell=True, check=False)
    print("BLAST done for", len(genomes), "genomes")


def load_hits(out: Path, g: str, query: str) -> pd.DataFrame:
    f = out / f"hits_all_{g}.txt"
    if f.stat().st_size == 0:
        return pd.DataFrame(columns=COLS)
    h = pd.read_csv(f, sep="\t", names=COLS)
    return h[h.qseqid.str.contains(query)]


def load_genome(data: Path, g: str) -> dict[str, str]:
    import gzip
    return {r.id: str(r.seq).upper()
            for r in SeqIO.parse(gzip.open(data / "genomes" / f"{g}.fna.gz", "rt"), "fasta")}


def clusters(hits: pd.DataFrame, cluster_gap: int = 500, evalue: float = 1e-5) -> list[dict]:
    """Group HSPs into (contig, proximity) clusters; return list sorted by summed bitscore."""
    hits = hits[hits.evalue < evalue].copy()
    if not len(hits):
        return []
    hits["smid"] = (hits[["sstart", "send"]].min(axis=1) + hits[["sstart", "send"]].max(axis=1)) / 2
    out = []
    for _, sub in hits.groupby("sseqid"):
        sub = sub.sort_values("smid")
        cur = [sub.iloc[0]]
        for _, row in sub.iloc[1:].iterrows():
            if abs(row.smid - cur[-1].smid) <= cluster_gap:
                cur.append(row)
            else:
                out.append(cur)
                cur = [row]
        out.append(cur)
    cl = []
    for c in out:
        cl.append(dict(contig=c[0].sseqid, bitscore=sum(r.bitscore for r in c),
                        smin=min(min(r.sstart, r.send) for r in c),
                        smax=max(max(r.sstart, r.send) for r in c),
                        qmin=min(r.qstart for r in c), qmax=max(r.qend for r in c),
                        minus=bool(c[0].sstart > c[0].send),
                        pident_best=max(r.pident for r in c), n_hsp=len(c)))
    return sorted(cl, key=lambda x: -x["bitscore"])


def extract_region(cl: dict, qlen: int, gd: dict[str, str], pad: int = 25) -> str:
    seq = gd[cl["contig"]]
    if cl["minus"]:
        lo = cl["smin"] - (qlen - cl["qmax"]) - pad
        hi = cl["smax"] + (cl["qmin"] - 1) + pad
    else:
        lo = cl["smin"] - (cl["qmin"] - 1) - pad
        hi = cl["smax"] + (qlen - cl["qmax"]) + pad
    lo, hi = max(lo, 1), min(hi, len(seq))
    sub = seq[lo - 1:hi]
    if cl["minus"]:
        sub = str(Seq(sub).reverse_complement())
    return sub


ALIGNER = PairwiseAligner()
ALIGNER.mode = "global"
ALIGNER.match_score = 2
ALIGNER.mismatch_score = -1
ALIGNER.open_gap_score = -2
ALIGNER.extend_gap_score = -0.5


def aln_metrics(query: str, subject: str):
    """Identity metrics with INTERNAL gaps only: padding columns at the subject
    extension are excluded by trimming to the span where the query aligns."""
    a = ALIGNER.align(query, subject)[0]
    s1, s2 = str(a[0]), str(a[1])
    qcols = [i for i, x in enumerate(s1) if x != "-"]
    lo, hi = qcols[0], qcols[-1]
    s1c, s2c = s1[lo:hi + 1], s2[lo:hi + 1]
    matches = sum(1 for x, y in zip(s1c, s2c) if x == y and x != "-" and y != "-")
    gaps = sum(1 for x, y in zip(s1c, s2c) if x == "-" or y == "-")
    return matches / len(query) * 100, matches / len(s1c) * 100, gaps, s1, s2


def main(data: str | Path, out: str | Path) -> None:
    data, out = Path(data), Path(out)
    out.mkdir(parents=True, exist_ok=True)

    ryhB, iscS_utr = load_mg1655_queries(data)
    (data / "queries.fasta").write_text(
        f">ryhB_MG1655\n{ryhB}\n>iscS_5UTR_MG1655\n{iscS_utr}\n")
    print("ryhB:", len(ryhB), "nt; iscS 5UTR+TIR:", len(iscS_utr), "nt; ATG check OK")

    genomes = sorted(Path(p).stem.replace(".fna", "")
                      for p in glob.glob(str(data / "genomes" / "*.fna.gz")))
    run_blast(data, out, genomes)

    rows, alignments = [], {}
    for g in genomes:
        gd = load_genome(data, g)
        entry = {"genome": g, "species": "E. coli" if g.startswith("E.coli") else "Salmonella"}
        rycl = clusters(load_hits(out, g, "ryhB"))
        ry_clusters = [c for c in rycl if c["qmax"] - c["qmin"] >= 30]
        entry["ryhB_n_loci"] = len(ry_clusters)
        aln_ry = []
        for i, c in enumerate(ry_clusters):
            sub = extract_region(c, 95, gd)
            idq, ida, gaps, _, _ = aln_metrics(ryhB, sub)
            aln_ry.append(dict(rank=i + 1, contig=c["contig"], smin=c["smin"], smax=c["smax"],
                                ident_qlen=round(idq, 1), ident_aln=round(ida, 1), gaps=gaps))
        entry["ryhB_loci"] = aln_ry
        entry["ryhB_ident"] = aln_ry[0]["ident_qlen"] if aln_ry else np.nan

        iscl = clusters(load_hits(out, g, "iscS"))
        iscl = [c for c in iscl if c["qmax"] - c["qmin"] >= 120]
        if iscl:
            sub = extract_region(iscl[0], 210, gd)
            idq, ida, gaps, s1, s2 = aln_metrics(iscS_utr, sub)
            entry.update(iscS_ident=round(idq, 1), iscS_ident_aln=round(ida, 1), iscS_gaps=gaps,
                         iscS_contig=iscl[0]["contig"], iscS_smin=iscl[0]["smin"], iscS_smax=iscl[0]["smax"])
            ry_sub = extract_region(ry_clusters[0], 95, gd) if aln_ry else None
            ry_s1, ry_s2 = (aln_metrics(ryhB, ry_sub)[3], aln_metrics(ryhB, ry_sub)[4]) if ry_sub else (None, None)
            alignments[g] = dict(iscS_q=s1, iscS_s=s2, ryhB_q=ry_s1, ryhB_s=ry_s2)
        else:
            entry.update(iscS_ident=np.nan, iscS_ident_aln=np.nan, iscS_gaps=np.nan)
        rows.append(entry)
        print(f'{g:26s} ryhB_loci={entry["ryhB_n_loci"]} best={entry["ryhB_ident"]}%  '
              f'iscS={entry.get("iscS_ident")}% (gaps={entry.get("iscS_gaps")})')

    flat = []
    for e in rows:
        best = e["ryhB_loci"][0] if e["ryhB_loci"] else {}
        flat.append(dict(genome=e["genome"], species=e["species"],
                          ryhB_n_loci=e["ryhB_n_loci"],
                          ryhB_ident=best.get("ident_qlen", np.nan),
                          ryhB_ident_aln=best.get("ident_aln", np.nan),
                          ryhB_gaps=best.get("gaps", np.nan),
                          ryhB_contig=best.get("contig"), ryhB_smin=best.get("smin"), ryhB_smax=best.get("smax"),
                          iscS_ident=e.get("iscS_ident"), iscS_ident_aln=e.get("iscS_ident_aln"),
                          iscS_gaps=e.get("iscS_gaps")))
    cons = pd.DataFrame(flat)
    cons.to_csv(out / "Table_S3_conservation_v3.csv", index=False)
    json.dump(rows, open(out / "conservation_paralogs.json", "w"), indent=1, cls=NpEncoder)
    json.dump(alignments, open(out / "alignments_v3.json", "w"), cls=NpEncoder)
    print("\n", cons[["genome", "ryhB_n_loci", "ryhB_ident", "ryhB_gaps", "iscS_ident", "iscS_gaps"]].to_string(index=False))

    L = len(iscS_utr)
    pos = np.full((len(genomes), L), np.nan)
    for gi, g in enumerate(genomes):
        a = alignments.get(g)
        if not a or not a.get("iscS_q"):
            continue
        qi = 0
        for x, y in zip(a["iscS_q"], a["iscS_s"]):
            if x != "-":
                pos[gi, qi] = 1.0 if (x == y and y != "-") else (0.0 if y != "-" else np.nan)
                qi += 1
    ecoli_idx = [i for i, g in enumerate(genomes) if g.startswith("E.coli")]
    salm_idx = [i for i, g in enumerate(genomes) if g.startswith("S.")]
    with np.errstate(invalid="ignore"):
        cons_all = np.nanmean(pos, axis=0) * 100
        cons_ec = np.nanmean(pos[ecoli_idx], axis=0) * 100
        cons_sal = np.nanmean(pos[salm_idx], axis=0) * 100
    coverage = np.mean(~np.isnan(pos), axis=0)

    win, half = 15, 7
    smooth = np.full(L, np.nan)
    for i in range(L):
        lo, hi = max(0, i - half), min(L, i + half + 1)
        if np.mean(~np.isnan(cons_all[lo:hi])) >= 0.7:
            smooth[i] = np.nanmean(cons_all[lo:hi])
    pd.DataFrame({"pos_rel_ATG": np.arange(L) - 150, "cons_all": cons_all, "cons_ecoli": cons_ec,
                  "cons_salmonella": cons_sal, "coverage": coverage, "smooth15": smooth}
                 ).to_csv(out / "iscS_perposition_conservation_v3.csv", index=False)

    rbs_lo, rbs_hi = 130, 149
    print(f"\nPairing/RBS window (ATG-20..-2): all={np.nanmean(cons_all[rbs_lo:rbs_hi]):.1f}%  "
          f"E.coli={np.nanmean(cons_ec[rbs_lo:rbs_hi]):.1f}%  Salmonella={np.nanmean(cons_sal[rbs_lo:rbs_hi]):.1f}%")
    print(f"Rest of UTR: all={np.nanmean(np.concatenate([cons_all[:rbs_lo], cons_all[rbs_hi:]])):.1f}%")
    print(f"First 60 nt coding: {np.nanmean(cons_all[150:]):.1f}%")

    dup = RNA.duplexfold(ryhB, iscS_utr)
    print(f"\nMG1655 reference RNAduplex MFE: {dup.energy:.1f} kcal/mol  structure: {dup.structure}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", required=True, help="directory with mg1655.fna.gz and genomes/*.fna.gz")
    p.add_argument("--out", required=True, help="output directory for BLAST hits and conservation tables")
    a = p.parse_args()
    main(a.data, a.out)
