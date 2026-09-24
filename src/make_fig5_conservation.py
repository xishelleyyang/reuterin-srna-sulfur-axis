#!/usr/bin/env python
"""Figure 5: conservation of ryhB and the iscS 5'UTR across the 13 target genomes.

Plotting logic is unchanged from the verified module-4 figure script (originally
developed in baseline/v3_working/make_figures_v3.py, lines ~526-687, then isolated
into its own script during the mSystems v4 finalization pass so there is exactly one
generating script per figure). Reads only the outputs of conservation_analysis.py.
"""
import argparse
import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def main(compgen_out, fig_out):
    os.makedirs(fig_out, exist_ok=True)
    OUT = compgen_out
    FIG = fig_out

    t3 = pd.read_csv(f"{OUT}/Table_S3_conservation_v3.csv")
    aln = json.load(open(f"{OUT}/alignments_v3.json"))
    aln = [dict(genome=k, **v) for k, v in aln.items()]
    perpos = pd.read_csv(f"{OUT}/iscS_perposition_conservation_v3.csv")


    def save(fig, name):
        fig.savefig(f"{FIG}/{name}.png", dpi=300, bbox_inches="tight")
        fig.savefig(f"{FIG}/{name}.svg", bbox_inches="tight")
        plt.close(fig)
        print("saved", name)


    def perpos_ident(qkey, skey, qlen):
        out = {}
        for e in aln:
            q, s = e[qkey], e[skey]
            ident = np.full(qlen, np.nan)
            qi = 0
            for qc, sc in zip(q, s):
                if qc != "-":
                    if qi < qlen:
                        ident[qi] = 1.0 if sc == qc else (0.0 if sc != "-" else np.nan)
                    qi += 1
            out[e["genome"]] = ident
        return out


    def smooth_edge(v, w=15):
        """Symmetric moving average with edge padding (no edge-decline artifact)."""
        v2 = pd.Series(v).interpolate(limit_direction="both").values
        k = np.ones(w) / w
        return np.convolve(np.pad(v2, w // 2, mode="edge"), k, mode="valid")


    EC = [e["genome"] for e in aln if e["genome"].startswith("E.coli")]
    SAL = [e["genome"] for e in aln if not e["genome"].startswith("E.coli")]

    ryh_pos = perpos_ident("ryhB_q", "ryhB_s", 95)
    ryh_ec = smooth_edge(np.nanmean([ryh_pos[g] for g in EC], axis=0) * 100)
    ryh_sal = smooth_edge(np.nanmean([ryh_pos[g] for g in SAL], axis=0) * 100)

    fig = plt.figure(figsize=(15.5, 10.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.34, wspace=0.24)

    ax = fig.add_subplot(gs[0, 0])
    df5 = t3.copy()
    df5["species"] = df5.genome.map(lambda g: "E. coli" if g.startswith("E.coli") else "Salmonella")
    df5 = df5.sort_values(["species", "ryhB_ident", "iscS_ident"]).reset_index(drop=True)
    names = [g.replace("E.coli_", "E. coli ").replace("S.", "S. ").replace("_", " ")
             for g in df5.genome]
    yy = np.arange(len(df5))
    for i, r in df5.iterrows():
        ax.plot([r.ryhB_ident, r.iscS_ident], [i, i], color="#CCCCCC", lw=2, zorder=1)
    ax.scatter(df5.ryhB_ident, yy, s=70, color="#D55E00", label="ryhB (95 nt)", zorder=3,
               edgecolor="k", linewidth=0.5)
    ax.scatter(df5.iscS_ident, yy, s=70, color="#0072B2", label="iscS 5'UTR+TIR (210 nt)",
               zorder=3, edgecolor="k", linewidth=0.5)
    ax.set_yticks(yy); ax.set_yticklabels(names, fontsize=9)
    for i, sp in enumerate(df5.species):
        if sp == "Salmonella":
            ax.get_yticklabels()[i].set_color("#7A4E00")
    ax.set_xlabel("Sequence identity to MG1655 (% of query length)", fontsize=11)
    ax.set_xlim(80, 102)
    ax.legend(fontsize=9, frameon=False, loc="lower left")
    ax.set_title("A  Conservation of ryhB and the iscS 5'UTR", loc="left", fontsize=12,
                 fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", alpha=0.3)
    ax.text(80.5, 5.6,
            "Salmonella: syntenic ryhB-1 shown; all 7 Salmonella\nalso carry a second, partial"
            " ryhB copy (core only)",
            fontsize=7.8, color="#7A4E00", va="top")

    axb = fig.add_subplot(gs[0, 1])
    xp = perpos.pos_rel_ATG.values
    cons_ec = smooth_edge(perpos.cons_ecoli.values)
    cons_sa = smooth_edge(perpos.cons_salmonella.values)
    axb.fill_between([-20, -2], 0, 105, color="#E9ED4C", alpha=0.5,
                     label="RyhB pairing site (iscS RBS)")
    axb.axvspan(0, 59, color="#EEEEEE", alpha=0.6, zorder=0)
    axb.plot(xp, cons_ec, color="#D55E00", lw=1.8, label="E. coli (n=6)")
    axb.plot(xp, cons_sa, color="#0072B2", lw=1.8, label="Salmonella (n=7)")
    for dpos in [-7, -74]:
        axb.annotate("", xy=(dpos, 99.2), xytext=(dpos, 103.5),
                     arrowprops=dict(arrowstyle="-|>", color="#0072B2", lw=1.4))
    axb.text(-22, 104.6, "Salmonella-shared\n1-nt deletions", fontsize=7.8, color="#0072B2",
             ha="center")
    axb.axvline(0, color="k", lw=1, ls="--")
    axb.text(2, 43.5, "iscS ATG", fontsize=9)
    axb.text(30, 103.5, "coding", fontsize=8, color="#666666", ha="center")
    axb.set_xlabel("Position relative to iscS start codon (nt)", fontsize=11)
    axb.set_ylabel("Identity to MG1655 (%)", fontsize=11)
    axb.set_ylim(40, 108)
    axb.set_xlim(-152, 62)
    axb.legend(fontsize=8.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16),
               ncol=3)
    axb.set_title("B  Position-resolved conservation of the iscS 5'UTR", loc="left",
                  fontsize=12, fontweight="bold")
    axb.spines[["top", "right"]].set_visible(False)

    axc = fig.add_subplot(gs[1, 0])
    axc.plot(np.arange(95), ryh_ec, color="#D55E00", lw=1.8, label="E. coli (n=6)")
    axc.plot(np.arange(95), ryh_sal, color="#0072B2", lw=1.8,
             label="Salmonella ryhB-1 (n=7)")
    axc.set_xlabel("Position along ryhB (nt)", fontsize=11)
    axc.set_ylabel("Identity to MG1655 (%)", fontsize=11)
    axc.set_ylim(60, 106)
    axc.legend(fontsize=8.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16),
               ncol=2)
    axc.set_title("C  Position-resolved conservation of ryhB", loc="left", fontsize=12,
                  fontweight="bold")
    axc.spines[["top", "right"]].set_visible(False)

    # (D) RBS alignment strip: query positions -20..-2 (0-based 130..148)
    # Salmonella rows on top so the shared-deletion annotation arrow stays clear
    axd = fig.add_subplot(gs[1, 1])
    lo_q, hi_q = 130, 149
    aln_d = ([e for e in aln if not e["genome"].startswith("E.coli")]
             + [e for e in aln if e["genome"].startswith("E.coli")])
    strip = []
    for e in aln_d:
        q, s = e["iscS_q"], e["iscS_s"]
        sub = []
        qi = 0
        for qc, sc in zip(q, s):
            if qc != "-":
                if lo_q <= qi < hi_q:
                    sub.append("." if sc == qc else (sc if sc != "-" else "-"))
                qi += 1
        strip.append(sub)
    strip = np.array(strip)
    qseq = []
    qi = 0
    for qc in aln_d[0]["iscS_q"]:
        if qc != "-":
            if lo_q <= qi < hi_q:
                qseq.append(qc)
            qi += 1
    genome_order = [e["genome"] for e in aln_d]
    labels5 = [g.replace("E.coli_", "E. coli ").replace("S.", "S. ").replace("_", " ")
               for g in genome_order]
    nrow = strip.shape[0]
    for i in range(nrow):
        for j in range(strip.shape[1]):
            ch = strip[i, j]
            fc = "#DFF0D8" if ch == "." else ("#4D4D4D" if ch == "-" else "#F4CCCC")
            axd.add_patch(plt.Rectangle((j, nrow - 1 - i), 1, 1, facecolor=fc,
                                        edgecolor="white", lw=0.5))
            if ch not in (".",):
                axd.text(j + 0.5, nrow - 1 - i + 0.5, ch, ha="center", va="center",
                         fontsize=7, color="k" if ch != "-" else "white", family="monospace")
    axd.set_xlim(0, strip.shape[1]); axd.set_ylim(0, nrow + 4.6)
    axd.set_xticks(np.arange(strip.shape[1]) + 0.5)
    axd.set_xticklabels(qseq, fontsize=7.5, family="monospace")
    axd.set_yticks(np.arange(nrow) + 0.5)
    axd.set_yticklabels(labels5[::-1], fontsize=8)
    for i, g in enumerate(genome_order[::-1]):
        if not g.startswith("E.coli"):
            axd.get_yticklabels()[i].set_color("#7A4E00")
    # annotate the Salmonella-shared deletion at query pos 143 (ATG -7) -> strip col 13
    del_col = 143 - lo_q
    axd.annotate("1-nt deletion shared by all 7 Salmonella (ATG −7)",
                 xy=(del_col + 0.5, nrow + 0.05), xytext=(del_col + 0.5, nrow + 2.4),
                 fontsize=8, color="#0072B2", ha="center",
                 arrowprops=dict(arrowstyle="-|>", color="#0072B2", lw=1.1))
    axd.set_xlabel("iscS RBS region (−20 to −2), MG1655 reference sequence below", fontsize=10)
    axd.set_title("D  RyhB pairing-site alignment across 13 genomes", loc="left", fontsize=12,
                  fontweight="bold")
    for sp in ["top", "right", "left", "bottom"]:
        axd.spines[sp].set_visible(False)
    axd.tick_params(length=0)
    leg = [Patch(facecolor="#DFF0D8", label="identical"), Patch(facecolor="#F4CCCC", label="mismatch"),
           Patch(facecolor="#4D4D4D", label="gap")]
    axd.legend(handles=leg, fontsize=8, frameon=False, loc="upper left",
               bbox_to_anchor=(0, -0.16), ncol=3)
    save(fig, "Fig5_conservation_v4")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--compgen-out", required=True, help="output directory from conservation_analysis.py")
    p.add_argument("--fig-out", required=True, help="directory to write Fig5_conservation_v4.{png,svg} (same directory make_figures.py uses for Fig1-4)")
    a = p.parse_args()
    main(a.compgen_out, a.fig_out)
    print("FIG5 DONE")
