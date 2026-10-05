#!/usr/bin/env python3
"""
make_plots.py — the figure set for a 750mer direct-print library QC run.

300 dpi PNG + matching SVG. Every figure is drawn from the merged metric bundle and
derived.json, so a figure cannot disagree with the report text.

Figures that carry a caveat draw it: the mispairing panel states the spec and its own axis
scale, the mosaic panel overlays measured detector sensitivity on the observed breakpoint
positions, and the deletion-rate comparison labels its denominator.
"""
import argparse, csv, json, os, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AP = argparse.ArgumentParser()
AP.add_argument("--metrics-dir", default=".")
AP.add_argument("--meta", required=True)
AP.add_argument("--outdir", default="plots")
AP.add_argument("--compare", default=None, help="optional prior-build metrics JSON for the head-to-head figure")
AP.add_argument("--microhomology", default=None)
A = AP.parse_args()
D = A.metrics_dir.rstrip("/") + "/"
F = A.outdir.rstrip("/") + "/"
os.makedirs(F, exist_ok=True)
J = lambda n: json.load(open(D + n))
_META = json.load(open(A.meta))
ORD = sorted({m["order"] for m in _META.values()})
SH = {o: o.split("-")[-1] for o in ORD}
DK, BL, OR, RD, GY, GR, PU = "#365F91", "#4F81BD", "#E8873A", "#C0392B", "#9AA5B1", "#4E8F5B", "#7B4B94"
_PAL = [RD, BL, GR, PU, OR, "#8A6D3B"]
OC = {o: _PAL[i % len(_PAL)] for i, o in enumerate(ORD)}
CLS = ["clean_designed", "within_well", "cross_well", "cross_sublib", "one_half", "ambiguous", "unmapped"]
CC = {"clean_designed": GR, "within_well": RD, "cross_well": OR, "cross_sublib": PU,
      "one_half": GY, "ambiguous": "#B0BFCE", "unmapped": "#D8D8D8"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.alpha": .25})
J = lambda n: json.load(open(D + n))


SKIPPED = []


def fig(name):
    """Decorator: draw a figure, or record why it was skipped. A missing optional input
    should cost one figure, not the whole run."""
    def deco(fn):
        def run(*a, **k):
            try:
                fn(*a, **k)
            except (KeyError, FileNotFoundError, IndexError, ValueError) as e:
                SKIPPED.append(f"{name}: {type(e).__name__} {e}")
                print(f"  SKIP {name} ({type(e).__name__}: {e})")
        run.__name__ = fn.__name__
        return run
    return deco


def save(fig, name):
    fig.tight_layout()
    fig.savefig(F + name + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(F + name + ".svg", bbox_inches="tight")
    plt.close(fig); print(" ", name)


LMIN = min(m["full_len"] for m in _META.values())
LMAX = max(m["full_len"] for m in _META.values())


def main():
    core, ext, der = J("metrics_core.json"), J("metrics_extra.json"), J("derived.json")
    BD, SC = J("blockdel_summary.json"), J("blockdel_scope.json")
    COV = {"cov_hist_10nt": SC.get("cov_hist_10nt_by_primer_status", {}),
           "cov_hist_both_by_sublib": SC.get("cov_hist_10nt_clean_by_sublib", {})}
    meta = _META

    # ---- fig01 ----
    try:
        fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
        for o in ORD:
            h = np.array(ext["rlen_hist_25bp_by_sublib_cls"].get(o + "|clean_designed", []), float)
            allh = collections.Counter()
            for c in ("clean_designed", "one_half", "unmapped"):
                for b, n in ext["rlen_hist_25bp_by_sublib_cls"].get(o + "|" + c, []): allh[b] += n
            xs = sorted(allh); ys = np.array([allh[x] for x in xs], float)
            ax[0].plot(xs, ys / ys.sum() * 100, color=OC[o], lw=1.4, label=f"{o} (clean+one-half+unmapped)")
            span_h = np.array(SC.get("mol_len_hist_10nt_gene_by_sublib", {}).get(o, []), float)
            if len(span_h):
                ax[1].plot(span_h[:, 0] + 5, span_h[:, 1] / span_h[:, 1].sum() * 100, color=OC[o], lw=1.4, label=o)
            else:
                ax[1].plot(h[:, 0], h[:, 1] / h[:, 1].sum() * 100, color=OC[o], lw=1.4, label=o)
        for a in ax:
            a.axvspan(LMIN, LMAX, color=GR, alpha=.12, zorder=0)
            a.set_ylabel("% of reads"); a.legend(fontsize=7, frameon=False)
        ax[0].set_xlabel("raw read length (nt, 25 nt bins)")
        ax[1].set_xlabel("primer-to-primer length (nt)")
        ax[0].set_xlim(0, 2000); ax[0].set_yscale("log")
        ax[0].set_title("All reads (log scale, raw length) — concatamers near 1,400 nt", loc="left", color=DK)
        ax[1].set_xlim(400, 900); ax[1].set_title(f"Both-primer primer-span — design {LMIN}–{LMAX} nt shaded", loc="left", color=DK)
        save(fig, "fig01_read_length")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig01' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig01:", e)

    # ---- fig02 ----
    try:
        fig, ax = plt.subplots(1, 2, figsize=(11.5, 3.8))
        x = np.arange(len(ORD)); bot = np.zeros(len(ORD))
        for c in CLS:
            v = np.array([core["per_sublib"][o].get(c + "_pct_of_all", 0) or 0 for o in ORD])
            ax[0].bar(x, v, bottom=bot, color=CC[c], label=c, width=.6)
            for i, (vv, bb) in enumerate(zip(v, bot)):
                if vv > 3: ax[0].text(i, bb + vv / 2, f"{vv:.1f}", ha="center", va="center", fontsize=7.5, color="white", fontweight="bold")
            bot += v
        ax[0].set_xticks(x); ax[0].set_xticklabels(ORD); ax[0].set_ylabel("% of ALL reads"); ax[0].set_ylim(0, 101)
        ax[0].legend(fontsize=7, ncol=2, loc="lower center", frameon=True)
        ax[0].set_title("Read composition (separate flowcells)", loc="left", color=DK)
        keys = ["within_well", "cross_well", "cross_sublib"]
        w = .25
        for j, k in enumerate(keys):
            v = [core["per_sublib"][o][k + "_pct_of_assigned"] for o in ORD]
            ax[1].bar(x + (j - 1) * w, v, w, color=CC[k], label=k)
            for i, vv in enumerate(v): ax[1].text(i + (j - 1) * w, vv * 1.03 + .001, f"{vv:.4f}", ha="center", fontsize=6.5, rotation=90)
        ax[1].set_xticks(x); ax[1].set_xticklabels(ORD); ax[1].set_ylabel("% of reads with BOTH halves assigned")
        ax[1].legend(fontsize=7.5, frameon=False)
        ax[1].set_ylim(0, max(core["per_sublib"][o][k + "_pct_of_assigned"] for o in ORD for k in keys) * 1.6)
        ax[1].set_title("Mispairing classes (spec ≤5% — note axis scale)", loc="left", color=DK)
        save(fig, "fig02_composition")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig02' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig02:", e)

    # ---- fig03 ----
    try:
        lc = ext["length_classes"]; rl = core["read_length_by_class"]
        cls = [c for c in CLS if c in rl]
        fig, ax = plt.subplots(figsize=(9.5, 3.6))
        for i, c in enumerate(cls):
            d = rl[c]
            ax.plot([d["p5"], d["p95"]], [i, i], color=CC[c], lw=3, solid_capstyle="butt")
            ax.plot(d["median"], i, "o", color=CC[c], ms=8)
            ax.text(max(d["p95"], 900) + 20, i, f"median {d['median']:.0f} · {lc[c]['pct_single']:.1f}% 600–900 nt · {lc[c]['pct_fused']:.1f}% ≥1,200 nt  (n={d['n']:,})", va="center", fontsize=7.5)
        ax.axvspan(LMIN, LMAX, color=GR, alpha=.14)
        ax.set_yticks(range(len(cls))); ax.set_yticklabels(cls); ax.invert_yaxis(); ax.set_xlim(0, 3300)
        ax.set_xlabel("read length (nt) · bar p5–p95, dot median")
        ax.set_title("Read length by class — a ligation-fused read sits near 2× the construct", loc="left", color=DK)
        save(fig, "fig03_readlen_by_class")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig03' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig03:", e)

    # ---- fig04 ----
    try:
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
        lab = {"both": "both primers", "fwd_only": "5' primer only", "tail_only": "3' primer only", "neither": "neither primer"}
        col = {"both": GR, "fwd_only": OR, "tail_only": BL, "neither": GY}
        for k in lab:
            h = np.array(COV["cov_hist_10nt"][k], float)
            if len(h): a1.plot(h[:, 0], h[:, 1] / h[:, 1].sum() * 100, color=col[k], lw=1.4, label=f"{lab[k]} (n={int(h[:, 1].sum()):,})")
        for o in ORD:
            h = np.array(COV["cov_hist_both_by_sublib"][o], float)
            a2.plot(h[:, 0], h[:, 1] / h[:, 1].sum() * 100, color=OC[o], lw=1.4, label=f"{o} ({SC['per_sublib_any_cls'][o]['pct']:.1f}% ≥50 nt deletion)")
        for a in (a1, a2):
            a.axvspan(LMIN, LMAX, color=GR, alpha=.14, zorder=0); a.set_yscale("log"); a.set_ylim(.005, 80); a.set_xlim(150, 780)
            a.set_xlabel("covered construct length (alignment span minus internal gaps, nt)"); a.legend(fontsize=7, frameon=False, loc="upper left")
        a1.set_ylabel("% of reads in class (log)")
        a1.set_title("By primer status, 1-in-8 subsample", loc="left", color=DK)
        a2.set_title("Both-primer reads, by sub-library", loc="left", color=DK)
        save(fig, "fig04_covered_length")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig04' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig04:", e)

    # ---- fig05 ----
    try:
        g = np.array(BD["gap_size_hist_10nt"], float)
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.3), gridspec_kw={"width_ratios": [1.3, 1]})
        a1.bar(g[:, 0], g[:, 1], width=9, color=RD); a1.axvline(50, color=DK, ls="--")
        a1.set_xlabel("deletion size (nt)"); a1.set_ylabel("reads (correctly paired, both primers)")
        a1.set_title("Block deletion size spectrum (≥50 nt calls)", loc="left", color=DK)
        ks = ["20", "30", "50", "100", "200"]
        a2.bar(range(5), [BD["block_deletion"][k]["pct"] for k in ks], color=BL, width=.6)
        for i, k in enumerate(ks): a2.text(i, BD["block_deletion"][k]["pct"] * 1.02, f"{BD['block_deletion'][k]['pct']:.2f}%", ha="center", fontsize=7.5)
        a2.set_xticks(range(5)); a2.set_xticklabels([f"≥{k}" for k in ks]); a2.set_xlabel("size threshold (nt)")
        a2.set_ylabel("% of both-primer clean reads"); a2.set_title("rate vs threshold (all orders)", loc="left", color=DK)
        save(fig, "fig05_deletion_spectrum")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig05' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig05:", e)

    # ---- fig06 ----
    try:
        ev = list(csv.DictReader(open(D + "blockdel_events.csv")))
        fig, ax = plt.subplots(1, 3, figsize=(12, 3.2), sharey=False)
        for a, o in zip(ax, ORD):
            e = [x for x in ev if x["sublib"] == o]
            st = [int(x["gap_start"]) for x in e]; en = [int(x["gap_start"]) + int(x["gap_len"]) for x in e]
            a.hist(st, bins=np.arange(0, 720, 10), color=RD, alpha=.7, label="start")
            a.hist(en, bins=np.arange(0, 720, 10), color=BL, alpha=.6, label="end")
            a.axvspan(175, 195, color=OR, alpha=.25, zorder=0)
            a.set_title(f"{o} (n={len(e):,} events)", loc="left", color=OC[o]); a.set_xlabel("position on construct (nt)")
            a.legend(fontsize=7, frameon=False)
        ax[0].set_ylabel("events"); ax[0].text(200, ax[0].get_ylim()[1] * .85, "internal forward-\nprimer site", fontsize=7, color=OR)
        save(fig, "fig06_deletion_position")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig06' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig06:", e)

    # ---- fig07 ----
    try:
        ds = der["design"]; att = der["blockdel"]["primer_site_attribution"]
        fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(12, 3.3))
        x = np.arange(3)
        a1.bar(x - .19, [ds[o]["pct_le3"] for o in ORD], .38, color=BL, label="internal site ≤3 mm (anywhere)")
        a1.bar(x + .19, [ds[o]["pct_site_175_195"] for o in ORD], .38, color=RD, label="site at nt 175–195")
        for i, o in enumerate(ORD): a1.text(i + .19, ds[o]["pct_site_175_195"] + 2, f"{ds[o]['pct_site_175_195']:.0f}%", ha="center", fontsize=8, color=RD)
        a1.set_xticks(x); a1.set_xticklabels([f"{SH[o]}\nFW {ds[o]['framework']}" for o in ORD]); a1.set_ylim(0, 100)
        a1.set_ylabel("% of designed constructs"); a1.legend(fontsize=6.8, frameon=False)
        a1.set_title("Design: forward-primer 3′ 12-mer internal match", loc="left", color=DK, fontsize=9)
        a2.bar(x, [att[o]["pct"] for o in ORD], color=OR, width=.55)
        for i, o in enumerate(ORD): a2.text(i, att[o]["pct"] + .5, f"{att[o]['pct']:.1f}%", ha="center", fontsize=8)
        a2.set_xticks(x); a2.set_xticklabels([SH[o] for o in ORD]); a2.set_ylabel("% of deletion events")
        a2.set_title("Observed: events start ≤40, end 170–200", loc="left", color=DK, fontsize=9)
        rates = [der["blockdel"]["per_gene_spread"][o] for o in ORD]
        a3.bar(x, [SC["per_sublib_any_cls"][o]["pct"] for o in ORD], color=[OC[o] for o in ORD], width=.55)
        for i, o in enumerate(ORD): a3.text(i, SC["per_sublib_any_cls"][o]["pct"] + .4, f"{SC['per_sublib_any_cls'][o]['pct']:.1f}%", ha="center", fontsize=8)
        a3.set_xticks(x); a3.set_xticklabels([SH[o] for o in ORD]); a3.set_ylabel("% both-primer reads with ≥50 nt deletion")
        a3.set_title("Deletion rate by order", loc="left", color=DK, fontsize=9)
        save(fig, "fig07_primer_collision")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig07' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig07:", e)

    # ---- fig08 ----
    try:
        pga = list(csv.DictReader(open(D + "blockdel_per_gene_all.csv")))
        fig, ax = plt.subplots(figsize=(8, 3.2))
        for o in ORD:
            r = np.array([int(z["g50"]) / int(z["n"]) * 100 for z in pga if z["sublib"] == o and int(z["n"]) >= 50])
            ax.hist(r, bins=np.arange(0, 101, 1.5), histtype="step", lw=1.6, color=OC[o], density=True,
                    label=f"{o}: median {np.median(r):.1f}% (n={len(r):,} variants ≥50 reads)")
        ax.set_xlabel("per-variant ≥50 nt deletion rate (%)"); ax.set_ylabel("density"); ax.set_xlim(0, 80)
        ax.legend(fontsize=7.5, frameon=False); ax.set_title("Deletion load per variant, by order", loc="left", color=DK)
        save(fig, "fig08_pervariant_deletion")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig08' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig08:", e)

    # ---- fig09 ----
    try:
        fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(12, 3.2))
        for k, c in (("fwd_only", OR), ("tail_only", BL)):
            h = np.array(BD["r_st_hist_10nt"][k], float); a1.plot(h[:, 0], h[:, 1] / h[:, 1].sum() * 100, color=c, label=lab[k])
            h = np.array(BD["r_en_hist_10nt"][k], float); a2.plot(h[:, 0], h[:, 1] / h[:, 1].sum() * 100, color=c, label=lab[k])
        a1.set_xlabel("alignment START on construct (nt)"); a2.set_xlabel("alignment END on construct (nt)")
        for a in (a1, a2): a.legend(fontsize=7, frameon=False); a.set_ylabel("% of class")
        a1.set_title("Where incomplete reads begin / end", loc="left", color=DK, fontsize=9)
        ks = ["both", "fwd_only", "tail_only", "neither"]
        fw = [100 * BD["strand"][k]["forward"] / max(1, BD["strand"][k]["forward"] + BD["strand"][k]["reverse"]) for k in ks]
        a3.bar(range(4), fw, color=GR, width=.55, label="forward"); a3.bar(range(4), [100 - v for v in fw], bottom=fw, color=GY, width=.55, label="reverse")
        for i, v in enumerate(fw): a3.text(i, v / 2, f"{v:.0f}%", ha="center", color="white", fontsize=8, fontweight="bold")
        a3.set_xticks(range(4)); a3.set_xticklabels([lab[k].replace(" primer", "") for k in ks], fontsize=7.5); a3.set_ylabel("% of reads")
        a3.legend(fontsize=7, frameon=False, loc="lower right"); a3.set_title("Strand by primer class", loc="left", color=DK, fontsize=9)
        save(fig, "fig09_truncation_strand")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig09' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig09:", e)

    # ---- fig10 ----
    try:
        fl = SC["fl_matrix"]; defs = ["aln_span_ge90", "readlen_ge90", "no_gap50", "readlen_and_nogap"]
        nm = {"aln_span_ge90": "alignment span ≥90%", "readlen_ge90": "primer-span length ≥90%", "no_gap50": "no internal deletion ≥50 nt",
              "readlen_and_nogap": "primer-span ≥90% AND no deletion"}
        dens = ["gene_assigned", "both_primer", "clean_designed"]
        dn = {"gene_assigned": "gene-assigned reads", "both_primer": "both-primer reads", "clean_designed": "correctly paired reads"}
        fig, ax = plt.subplots(figsize=(9, 3.4)); w = .26
        for j, dk in enumerate(dens):
            v = [fl[d][dk] for d in defs]
            ax.barh(np.arange(4) + (j - 1) * w, v, w, color=[GY, BL, GR][j], label=dn[dk])
            for i, vv in enumerate(v): ax.text(vv + .5, i + (j - 1) * w, f"{vv:.1f}", va="center", fontsize=6.8)
        ax.axvline(60, color=DK, ls="--"); ax.set_yticks(range(4)); ax.set_yticklabels([nm[d] for d in defs]); ax.invert_yaxis()
        ax.set_xlim(0, 112); ax.set_xlabel("% full-length"); ax.legend(fontsize=7, frameon=False, loc="lower left")
        ax.set_title("Full-length rate: the definition matters more than the denominator", loc="left", color=DK)
        save(fig, "fig10_full_length_definitions")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig10' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig10:", e)

    # ---- fig11 ----
    try:
        pg = {r["gene"]: int(r["n"]) for r in csv.DictReader(open(D + "per_gene.csv"))}
        U = der["uniformity"]
        fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.5))
        for o in ORD:
            n = np.array(sorted([pg.get(g, 0) for g in meta if meta[g]["order"] == o], reverse=True), float)
            ax[0].plot(np.arange(1, len(n) + 1) / len(n) * 100, n, color=OC[o], lw=1.2, label=f"{o}: median {U[o]['median']:.0f}, 95/5 {U[o]['p95p5']:.2f}")
            srt = np.sort(n); cum = np.cumsum(srt) / srt.sum()
            ax[1].plot(np.linspace(0, 1, len(cum)), cum, color=OC[o], lw=1.4, label=f"{o}: Gini {U[o]['gini']:.3f}, CV {U[o]['cv']:.3f}")
        ax[0].set_yscale("log"); ax[0].set_xlabel("variant rank (percentile)"); ax[0].set_ylabel("correctly paired reads"); ax[0].legend(fontsize=6.8, frameon=False)
        ax[1].plot([0, 1], [0, 1], ls=":", color=GY); ax[1].set_xlabel("cumulative share of variants"); ax[1].set_ylabel("cumulative share of reads"); ax[1].legend(fontsize=6.8, frameon=False)
        ax[0].set_title("Rank abundance", loc="left", color=DK); ax[1].set_title("Lorenz curve", loc="left", color=DK)
        rar = core.get("rarefaction") or []
        if not rar:
            raise KeyError("rarefaction (depth subsampling not run)")
        ax[2].plot([r["n_reads"] for r in rar], [r["n_variants"] for r in rar], "o-", color=DK)
        ax[2].axhline(core["per_variant"]["designed"], color=OR, ls="--"); ax[2].set_xlabel("correctly paired reads sampled (all orders)")
        ax[2].set_ylabel("distinct variants seen"); ax[2].set_title("Rarefaction", loc="left", color=DK)
        ax[2].set_ylim(core["per_variant"]["designed"] * .98, core["per_variant"]["designed"] * 1.003)
        save(fig, "fig11_uniformity")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig11' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig11:", e)

    # ---- fig12 ----
    try:
        ROWS = "ABCDEFGHIJKLMNOP"
        pw = list(csv.DictReader(open(D + "per_well.csv")))
        pwd = der["blockdel"]["per_well"]
        fig, ax = plt.subplots(3, 3, figsize=(13, 9.5))
        panels = [("depth", "correctly paired reads / well", "viridis"), ("within", "within-well chimera %", "magma_r"), ("del", "≥50 nt deletion %", "Reds")]
        for ci, o in enumerate(ORD):
            for ri, (k, lbl, cmap) in enumerate(panels):
                g = np.full((16, 12), np.nan)
                for r in pw:
                    if r["order"] != o: continue
                    w = r["well"]; y, xx = ROWS.index(w[0]), int(w[1:]) - 1
                    if k == "depth": g[y, xx] = int(r["n_clean"])
                    elif k == "within": g[y, xx] = 100 * int(r["n_within"]) / int(r["n_assigned"])
                    else: g[y, xx] = pwd[o]["_wells"].get(w, np.nan)
                a = ax[ri, ci]
                vm = {"within": (0, .1), "del": (0, max(max(pwd[z]["_wells"].values()) for z in ORD))}.get(k, (None, None))
                im = a.imshow(g, cmap=cmap, aspect="auto", vmin=vm[0], vmax=vm[1])
                a.set_xticks(range(12)); a.set_xticklabels(range(1, 13), fontsize=6); a.set_yticks(range(16)); a.set_yticklabels(list(ROWS), fontsize=6)
                a.grid(False); plt.colorbar(im, ax=a, fraction=.046).set_label(lbl, fontsize=7)
                if ri == 0: a.set_title(o, color=OC[o], fontweight="bold")
        fig.suptitle("Plate maps (replicate-1 positions, cols 1–12; replicate 2 in cols 13–24 is sequence-identical and not separable)", fontsize=10, color=DK)
        save(fig, "fig12_plate_maps")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig12' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig12:", e)

    # ---- fig13 ----
    try:
        rec = der["recoverability"]
        fig, ax = plt.subplots(figsize=(8.5, 3.4)); x = np.arange(3); w = .13
        for j, t in enumerate(("90", "95", "99")):
            ax.bar(x + (j - 1) * 2 * w - w / 2, [rec[o][t]["ignoring_deletions"] for o in ORD], w, color=GY, label="ignoring deletions" if j == 0 else None)
            ax.bar(x + (j - 1) * 2 * w + w / 2, [rec[o][t]["counting_deletions"] for o in ORD], w, color=[OC[o] for o in ORD], label=None)
            for i, o in enumerate(ORD): ax.text(i + (j - 1) * 2 * w + w / 2, rec[o][t]["counting_deletions"] * 1.03, f"{t}%", ha="center", fontsize=6.5)
        ax.set_xticks(x); ax.set_xticklabels(ORD); ax.set_ylabel("cells screened (minimum)"); ax.legend(fontsize=7.5, frameon=False)
        ax.set_title("Cells to recover 90 / 95 / 99% of variants as intact, correctly paired clones (per sub-library)", loc="left", color=DK, fontsize=9)
        save(fig, "fig13_recoverability")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig13' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig13:", e)

    # ---- fig14 ----
    try:
        if not A.compare:
            raise FileNotFoundError("no --compare metrics supplied")
        a13 = json.load(open(A.compare))
        c13 = a13["metrics_core"]["per_sublib"]; b13 = a13["blockdel_per_sublib"]
        s13 = a13.get("blockdel_scope", {})
        CMP = sorted(c13)
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.4))
        labs = [o.split("-")[-1] + "\n" + a13.get("frameworks", {}).get(o, "") for o in CMP] \
             + [SH[o] + "\n" + (der.get("design", {}).get(o, {}).get("framework", "")) for o in ORD]
        wv = [c13[o]["within_well_pct_of_assigned"] for o in CMP] + [core["per_sublib"][o]["within_well_pct_of_assigned"] for o in ORD]
        dv = [b13[o]["pct"] for o in CMP] + [BD["per_sublib"][o]["pct"] for o in ORD]
        n = len(wv)
        cols = [GY] * len(CMP) + [OC[o] for o in ORD]
        a1.bar(range(n), wv, color=cols); a2.bar(range(n), dv, color=cols)
        for i in range(n):
            a1.text(i, wv[i] * 1.02, f"{wv[i]:.4f}", ha="center", fontsize=7)
            a2.text(i, dv[i] * 1.02, f"{dv[i]:.1f}", ha="center", fontsize=7)
        for a in (a1, a2): a.set_xticks(range(n)); a.set_xticklabels(labs, fontsize=7.5)
        a1.set_ylabel("% of assigned reads"); a1.set_title("Within-well pairing chimera", loc="left", color=DK)
        a2.set_ylabel("% of both-primer correctly paired reads"); a2.set_title("≥50 nt internal block deletion", loc="left", color=DK)
        fig.text(.5, -.02, "Deletion rates on both panels: correctly paired, both-primer reads, 1-in-8 subsample, identical pipeline for both libraries.", ha="center", fontsize=7, color=GY)
        save(fig, "fig14_vs_ailk013")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig14' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig14:", e)

    # ---- fig15 ----
    try:
        if os.path.exists(D + "mosaic_raw.json"):
            mos = J("mosaic_raw.json"); ctl = J("controls.json")["mospos"]
            fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.3))
            h = mos["position_hist"]; ctr = np.arange(.05, 1, .1)
            a1.bar(ctr, h, width=.09, color=BL, label="observed breakpoints")
            b = a1.twinx(); sx = [.1, .2, .3, .4, .5, .6, .7, .8, .9]
            b.plot(sx, [ctl[f"flat_{p:g}"]["pct"] for p in sx], "o-", color=OR); b.set_ylim(0, 105); b.set_ylabel("detector sensitivity %", color=OR); b.grid(False)
            a1.set_xlabel("breakpoint position (fraction of read)"); a1.set_ylabel("breakpoints")
            a1.set_title(f"Intra-segment breakpoints ({sum(h)} events) vs sensitivity", loc="left", color=DK, fontsize=9)
            vals = [mos["mosaic_pct_with_margin"], mos["mosaic_pct_no_margin"]]
            a2.bar(["vote margin 3\n(shipped)", "no margin\n(control)"], vals, color=[BL, GY], width=.5)
            for i, v in enumerate(vals): a2.text(i, v * 1.03, f"{v:.3f}%", ha="center")
            a2.set_ylabel("% of full-length reads with ≥1 breakpoint"); a2.set_title(f"Mosaic rate (n={mos['n_full_length_sampled']:,} full-length reads)", loc="left", color=DK, fontsize=9)
            save(fig, "fig15_mosaic")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig15' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig15:", e)

    # ---- fig16 ----
    try:
        if not A.microhomology:
            raise FileNotFoundError("no --microhomology table supplied")
        mh = {(r["order"], r["well"]): r for r in csv.DictReader(open(A.microhomology))}
        fig, ax = plt.subplots(1, 3, figsize=(12, 3.2))
        for a, o in zip(ax, ORD):
            xs, ys = [], []
            for r in pw:
                k = (r["order"], r["well"])
                if r["order"] != o or k not in mh or int(mh[k]["n_mates"]) != 120: continue
                xs.append(int(mh[k]["shared_18mers"])); ys.append(100 * int(r["n_within"]) / int(r["n_assigned"]))
            rr = der["regression"][o]["shared_18mers"]
            a.scatter(xs, ys, s=8, alpha=.5, color=OC[o]); a.set_xlabel("shared 18-mers in well (design)")
            a.set_title(f"{o}: r = {rr['r']:+.3f} ({rr['ci'][0]:+.2f}, {rr['ci'][1]:+.2f}), n={rr['n']}", loc="left", fontsize=8.5, color=OC[o])
        ax[0].set_ylabel("within-well chimera % (per well)")
        save(fig, "fig16_microhomology")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append('fig16' + ": " + type(e).__name__ + " " + str(e)); print("  SKIP fig16:", e)

    # ---- fig17 ----
    try:
        # 17a / customer collateral: primer-to-primer span, not raw FASTQ length.
        # BLOCKDEL is a 1-in-N subsample; scale counts up so the axis stays in millions.
        mol = SC.get("mol_len_hist_10nt") or []
        if not mol:
            raise KeyError("mol_len_hist_10nt")
        h = np.array(mol, float)
        sub = int(SC.get("subsample") or 8)
        ys = h[:, 1] * sub / 1e6
        def readlen_span(ax, title=True):
            ax.bar(h[:, 0] + 5, ys, width=9, color=BL)
            ax.axvspan(LMIN, LMAX, color=GR, alpha=.18, zorder=0)
            top = float(ys.max()) if len(ys) else 1.0
            ax.annotate(f"designed window\n{LMIN}–{LMAX} bp", (LMAX, top * .95), xytext=(1050, top * .85), color=DK, fontsize=8,
                        arrowprops=dict(arrowstyle="-", color=DK, lw=.9))
            ax.set_xlim(0, 2000); ax.set_xlabel("Primer-to-primer length (bp)"); ax.set_ylabel("Reads (millions)")
            if title:
                n_span = (SC.get("primer_span") or {}).get("n_span") or int(h[:, 1].sum())
                ax.set_title(f"Primer-to-primer length, both-primer reads (n≈{n_span * sub:,}; 10 nt bins; ONT adapters excluded)", loc="left", color=DK)
        fig, ax = plt.subplots(figsize=(8, 3.2)); readlen_span(ax); save(fig, "fig17a_readlen_all_linear")
        fig, ax = plt.subplots(figsize=(6.2, 2.4)); readlen_span(ax, title=False); save(fig, "collateral_readlen")

        # 17b: per sub-library, linear, 0-2000, stacked by class
        cols = {"clean_designed": GR, "one_half": "#9AA5B1", "unmapped": "#D8D8D8"}
        fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.3), sharey=True)
        for a, o in zip(axs, ORD):
            bot = None
            for c in ("clean_designed", "one_half", "unmapped"):
                hh = dict((int(b), n) for b, n in ext["rlen_hist_25bp_by_sublib_cls"].get(o + "|" + c, []))
                xs = np.arange(0, 2025, 25); ys = np.array([hh.get(int(x), 0) for x in xs], float) / 1e6
                a.bar(xs + 12.5, ys, width=24, bottom=bot, color=cols[c], label=c)
                bot = ys if bot is None else bot + ys
            a.axvspan(LMIN, LMAX, color=GR, alpha=.15, zorder=0); a.set_xlim(0, 2000)
            a.set_title(f"{o} ({core['per_sublib'][o]['n_reads']/1e6:.1f} M reads)", loc="left", color=OC[o]); a.set_xlabel("Read length (bp)")
        axs[0].set_ylabel("Reads (millions)"); axs[0].legend(fontsize=7, frameon=False)
        fig.suptitle("Read length by sub-library and class (within-well / cross-well / cross-sub-library reads are <0.5% and not drawn)", color=DK, fontsize=9.5)
        save(fig, "fig17b_readlen_by_sublib")

        # 17c: zoom 400-900, primer-span vs designed-length distribution
        fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.2), sharey=False)
        by = SC.get("mol_len_hist_10nt_gene_by_sublib", {})
        pspan = SC.get("primer_span_by_sublib", {})
        for a, o in zip(axs, ORD):
            L = [m["full_len"] for m in _META.values() if m["order"] == o]
            if o in by and by[o]:
                hh = np.array(by[o], float)
                sel = (hh[:, 0] >= 400) & (hh[:, 0] < 900)
                tot = hh[:, 1].sum()
                a.bar(hh[sel, 0] + 5, hh[sel, 1] / tot * 100, width=9, color=OC[o], alpha=.75,
                      label="primer-span (both primers)")
                c = collections.Counter(int(x // 10 * 10) for x in L)
                xs = np.arange(400, 900, 10)
                a.plot(xs + 5, [c.get(int(x), 0) / len(L) * 100 for x in xs], "o-", color="k", lw=1.1, ms=3,
                       label="designed constructs")
            else:
                hh = np.array(ext["rlen_hist_25bp_by_sublib_cls"][o + "|clean_designed"], float)
                sel = (hh[:, 0] >= 400) & (hh[:, 0] < 900)
                a.bar(hh[sel, 0] + 12.5, hh[sel, 1] / hh[:, 1].sum() * 100, width=24, color=OC[o], alpha=.75,
                      label="raw read length")
                c = collections.Counter(int(x // 25 * 25) for x in L)
                xs = np.arange(400, 900, 25)
                a.plot(xs + 12.5, [c.get(int(x), 0) / len(L) * 100 for x in xs], "o-", color="k", lw=1.1, ms=3,
                       label="designed constructs")
            med = (pspan.get(o) or {}).get("median")
            a.set_title(f"{o}: primer-span median {med} nt (design median {int(np.median(L))} nt)" if med
                        else f"{o}: design median {int(np.median(L))} nt", loc="left", fontsize=8.5, color=OC[o])
            a.set_xlabel("length (bp)"); a.set_ylabel("% of reads / % of designs"); a.set_xlim(400, 900)
            a.legend(fontsize=7, frameon=False, loc="upper left")
        fig.suptitle("Zoom 400–900 bp: primer-to-primer span vs designed construct length (ONT adapters outside the primers excluded)", color=DK, fontsize=9.5)
        save(fig, "fig17c_readlen_zoom_vs_design")
    except (KeyError, FileNotFoundError, IndexError, ValueError, TypeError) as e:
        SKIPPED.append("fig17: " + type(e).__name__ + " " + str(e)); print("  SKIP fig17:", e)

    if SKIPPED:
        json.dump(SKIPPED, open(F + "figures_skipped.json", "w"), indent=1)
        print("skipped figures:", len(SKIPPED))


if __name__ == "__main__":
    main()
