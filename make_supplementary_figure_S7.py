import os

import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr

# ─── USER SETTINGS ────────────────────────────────────────────────────────────
XLSX_PATH = (
    "data/Excel/ml_assessment_library/EVAAV/enrichment_score_with_amino_acid_sequence.xlsx"
)
PLOTS_DIR = "plots"

REP_COLS = ["enrichment_score_rep_1", "enrichment_score_rep_2", "enrichment_score_rep_3"]

# ML-designed sequences to look up (variants of the scaffold with Cys walked
# across positions 2-7, plus an all-Cys and a double-Cys design)
QUERY_SEQUENCES = [
    "GCCCCCCS",
    "GCGGSGGS",
    "GGCGSGGS",
    "GGSCSGGS",
    "GGSGCGGS",
    "GGSGGCGS",
    "GGSGGSCS",
    "GCGGSGCG",
]

SEQ_LEN = 8  # all query sequences are 8 residues long; restrict library comparisons to this length

# Wild-type loop sequence (pre-randomization), per BenchmarkingMLvsDEtop10 - Reviews.ipynb
WT_SEQ = "GSGQNQQ"

# ─── STYLE ────────────────────────────────────────────────────────────────────
BASE_FS = 16
rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "axes.titlesize": BASE_FS,
    "axes.labelsize": BASE_FS,
    "xtick.labelsize": BASE_FS - 2,
    "ytick.labelsize": BASE_FS - 2,
    "legend.fontsize": BASE_FS - 2,
})

COL_NO_CYS = "#999999"    # grey  - no cysteine
COL_CYS = "#CC3311"       # red   - has cysteine
COL_QUERY = "#332288"     # indigo - ML-designed query sequences
COL_MEDIAN = "#333333"    # near-black reference line
COL_WT = "#117733"        # green - wild-type baseline


# ─── WILD-TYPE BASELINE ───────────────────────────────────────────────────────
def get_wt_baseline(df):
    """Wild-type loop sequence's enrichment score, averaged over all its DNA variants."""
    matches = df[df["amino_acid_sequence"] == WT_SEQ]
    if matches.empty:
        raise ValueError(f"Wild-type sequence {WT_SEQ!r} not found in {XLSX_PATH}")
    mean = matches["average_enrichment_score"].mean()
    sem = matches["average_enrichment_score"].sem()
    n = len(matches)
    print(
        f"Wild-type baseline {WT_SEQ!r} (n_Cys={WT_SEQ.count('C')}): "
        f"average_enrichment_score = {mean:.3f} +/- {sem:.3f} SEM (n={n} DNA variants)"
    )
    return {"sequence": WT_SEQ, "mean": mean, "sem": sem, "n": n}


# ─── LOOKUP ───────────────────────────────────────────────────────────────────
def lookup_sequences(df, sequences):
    """Look up each amino acid sequence and report its enrichment score(s)."""
    rows = []
    print("=" * 78)
    print("Sequence lookup")
    print("=" * 78)
    df_len = df[df["amino_acid_sequence"].str.len() == SEQ_LEN]
    for seq in sequences:
        matches = df[df["amino_acid_sequence"] == seq]
        n_cys = seq.count("C")
        if matches.empty:
            print(f"  {seq}  (n_Cys={n_cys})  -> NOT FOUND in library")
            rows.append({
                "amino_acid_sequence": seq, "n_cys": n_cys, "found": False,
                "n_dna_variants": 0, "average_enrichment_score": np.nan,
                "std_dev_enrichment_score": np.nan, "percentile_in_library": np.nan,
            })
            continue

        avg = matches["average_enrichment_score"].mean()
        std = matches["std_dev_enrichment_score"].mean()
        percentile = (df_len["average_enrichment_score"] < avg).mean() * 100
        print(
            f"  {seq}  (n_Cys={n_cys})  -> average_enrichment_score = {avg:.3f} "
            f"+/- {std:.3f}  (n={len(matches)} DNA variant(s), "
            f"{percentile:.1f}th percentile of {SEQ_LEN}-mer library)"
        )
        rows.append({
            "amino_acid_sequence": seq, "n_cys": n_cys, "found": True,
            "n_dna_variants": len(matches),
            "average_enrichment_score": avg,
            "std_dev_enrichment_score": std,
            "percentile_in_library": percentile,
        })
    print("=" * 78)
    return pd.DataFrame(rows)


# ─── LIBRARY-WIDE CYS ANALYSIS ────────────────────────────────────────────────
def build_library_table(df):
    """One row per unique 8-residue amino acid sequence, averaged across DNA variants."""
    df_len = df[df["amino_acid_sequence"].str.len() == SEQ_LEN].copy()
    lib = (
        df_len.groupby("amino_acid_sequence", as_index=False)["average_enrichment_score"]
        .mean()
    )
    lib["n_cys"] = lib["amino_acid_sequence"].str.count("C")
    for pos in range(SEQ_LEN):
        lib[f"pos_{pos + 1}_is_cys"] = lib["amino_acid_sequence"].str[pos] == "C"
    return lib


def plot_enrichment_by_cys_count(lib, lookup_df, wt, out_path):
    """Fig 1: distribution of enrichment score across the library, grouped by Cys count,
    shown as box-and-whisker plots, with the ML-designed query sequences overlaid."""
    bins = [0, 1, 2, 3]
    labels = ["0", "1", "2", "3+"]
    lib = lib.copy()
    lib["cys_bin"] = pd.cut(
        lib["n_cys"], bins=[-0.5, 0.5, 1.5, 2.5, 100], labels=labels
    )

    # extra spacing between box centers leaves room for the "med=" labels
    # beside each box without overlapping the neighboring box
    SPACING = 1.3

    fig, ax = plt.subplots(figsize=(8, 6))
    box_data, box_bins, box_colors = [], [], []
    for i, label in enumerate(labels):
        vals = lib.loc[lib["cys_bin"] == label, "average_enrichment_score"]
        if len(vals) == 0:
            continue
        box_data.append(vals)
        box_bins.append(i)
        box_colors.append(COL_NO_CYS if label == "0" else COL_CYS)
    box_positions = [b * SPACING for b in box_bins]

    bp = ax.boxplot(
        box_data, positions=box_positions, widths=0.5, patch_artist=True,
        showfliers=False,
        medianprops=dict(color=COL_MEDIAN, lw=2.5),
        whiskerprops=dict(color="0.3"), capprops=dict(color="0.3"),
        boxprops=dict(edgecolor="0.3"),
        zorder=2,
    )
    for patch, color in zip(bp["boxes"], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    for pos, vals in zip(box_positions, box_data):
        median = vals.median()
        ax.text(pos + 0.32, median, f"med={median:.2f}", va="center", fontsize=BASE_FS - 4)

    # overlay the found ML-designed query sequences, numbered to avoid overlapping labels
    found = lookup_df[lookup_df["found"]].reset_index(drop=True)
    query_bins = found["n_cys"].clip(upper=3).map({0: 0, 1: 1, 2: 2, 3: 3})
    x_jitter_q = query_bins.astype(float) * SPACING
    for b in query_bins.unique():
        idx = query_bins[query_bins == b].index
        if len(idx) > 1:
            offsets = np.linspace(-0.22, 0.22, len(idx))
            for offset, i in zip(offsets, found.loc[idx].sort_values("average_enrichment_score").index):
                x_jitter_q[i] = b * SPACING + offset
    ax.scatter(
        x_jitter_q, found["average_enrichment_score"], s=110, marker="D",
        color=COL_QUERY, edgecolor="white", linewidths=1.2, zorder=5,
        label="Cysteine-linker sequence controls",
    )
    for i, (_, row) in enumerate(found.iterrows()):
        ax.annotate(
            str(i + 1), (x_jitter_q[i], row["average_enrichment_score"]),
            ha="center", va="center", fontsize=BASE_FS - 7, color="white", fontweight="bold", zorder=6,
        )
    legend_text = "\n".join(f"{i + 1}: {s}" for i, s in enumerate(found["amino_acid_sequence"]))
    ax.text(
        0.98, 0.98, legend_text, transform=ax.transAxes, ha="right", va="top",
        fontsize=BASE_FS - 5, color=COL_QUERY, family="monospace",
        bbox=dict(facecolor="white", edgecolor=COL_QUERY, linewidth=0.8, alpha=0.9, boxstyle="round,pad=0.4"),
    )

    ax.axhline(wt["mean"], color=COL_WT, linestyle="--", lw=1.8, zorder=3,
               label=f"Wild-type ({wt['sequence']}, {wt['mean']:.2f})")

    # y-max: highest non-outlier whisker top across the boxes (plus query
    # sequences / WT line), with 15% headroom, rounded up to the nearest 0.5
    def _whisker_top(vals, whis=1.5):
        q1, q3 = vals.quantile(0.25), vals.quantile(0.75)
        within_fence = vals[vals <= q3 + whis * (q3 - q1)]
        return within_fence.max() if len(within_fence) else vals.max()

    y_max_data = max(
        [_whisker_top(vals) for vals in box_data]
        + [found["average_enrichment_score"].max(), wt["mean"]]
    )
    y_max = np.ceil(y_max_data * 1.15 * 2) / 2

    ax.set_xlabel("Number of cysteine residues")
    ax.set_xticks([b * SPACING for b in range(len(labels))])
    ax.set_xticklabels(labels)
    ax.set_xlim(-0.6, (len(labels) - 1) * SPACING + 0.9)
    ax.set_ylim(0, y_max)
    ax.set_ylabel(
        r"Average EV enrichment score $S = \frac{f_{\mathrm{EV\!-\!AAV}}}{f_{\mathrm{plasmid}}}$"
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"Saved: {out_path}")

def main():
    os.makedirs(PLOTS_DIR, exist_ok=True)
    df = pd.read_excel(XLSX_PATH)

    lookup_df = lookup_sequences(df, QUERY_SEQUENCES)

    wt = get_wt_baseline(df)
    lib = build_library_table(df)

    rho, rho_p = spearmanr(lib["n_cys"], lib["average_enrichment_score"])
    print(
        f"\nLibrary-wide correlation (n_Cys vs. average_enrichment_score, "
        f"n={len(lib)} unique {SEQ_LEN}-mers): Spearman rho = {rho:.3f}, p = {rho_p:.2e}"
    )

    plot_enrichment_by_cys_count(lib, lookup_df, wt, os.path.join(PLOTS_DIR, "supplementary_figure_S7.png"))

if __name__ == "__main__":
    main()
