"""Compares per-subfamily TE counts of bowtie2 alignment modes against -k 100

Usage: python compare_modes.py <analysis directory>
Writes te_alignment_comparison/comparison/*.tsv and prints a report.
"""

import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

C = os.path.join(sys.argv[1], "te_alignment_comparison")
MODES = ["default", "very_sensitive", "k10", "k20", "k100"]
REF = "k100"
# Young L1 and IAP subfamilies
YOUNG = [
    "L1MdTf_I",
    "L1MdTf_II",
    "L1MdTf_III",
    "L1MdA_I",
    "L1MdA_II",
    "L1MdA_III",
    "L1MdGf_I",
    "L1MdGf_II",
    "L1MdN_I",
    "L1MdF_I",
    "L1MdF_II",
    "IAPEz-int",
    "IAPLTR1_Mm",
    "IAPLTR1a_Mm",
]

out_dir = os.path.join(C, "comparison")
os.makedirs(out_dir, exist_ok=True)
modes = [m for m in MODES if os.path.exists(f"{C}/counts/{m}.multi.txt")]


def read_counts(mode, kind):
    df = pd.read_csv(f"{C}/counts/{mode}.{kind}.txt", sep="\t", comment="#")
    df = df.set_index("Geneid").iloc[:, 5:]
    df.columns = [os.path.basename(c).split(".")[0] for c in df.columns]
    return df


def read_time(sample, mode):
    with open(f"{C}/logs/time_{sample}_{mode}.txt") as f:
        m = re.search(r"([\d.]+) s wall, (\d+) kB", f.read())
    return float(m.group(1)) / 60, int(m.group(2)) / 1024**2


def cpm(df):
    return df / df.sum() * 1e6


multi = {m: read_counts(m, "multi") for m in modes}
unique = {m: read_counts(m, "unique") for m in modes}
ref = multi[REF]
samples = list(ref.columns)

# Run time and peak memory
rows = []
for m in modes:
    for s in samples:
        minutes, gb = read_time(s, m)
        rows.append(
            {
                "mode": m,
                "sample": s,
                "minutes": round(minutes, 1),
                "max_GB": round(gb, 1),
            }
        )
times = pd.DataFrame(rows)
times.to_csv(f"{out_dir}/run_time.tsv", sep="\t", index=False)
print("## Run time and peak memory")
print(times.pivot(index="mode", columns="sample", values="minutes").loc[modes])
print(times.groupby("mode")["max_GB"].max().loc[modes])

# Agreement with the reference mode (subfamilies with >= 50 counts)
rows = []
for m in modes:
    if m == REF:
        continue
    for s in samples:
        keep = ref[s] >= 50
        x, y = multi[m].loc[keep, s], ref.loc[keep, s]
        lr = np.log2(((x + 1) / x.sum()) / ((y + 1) / y.sum()))
        rows.append(
            {
                "mode": m,
                "sample": s,
                "n": keep.sum(),
                "pearson_log2": round(pearsonr(np.log2(x + 1), np.log2(y + 1))[0], 4),
                "spearman": round(spearmanr(x, y)[0], 4),
                "median_abs_log2": round(lr.abs().median(), 3),
                "n_abs_log2_gt_0.5": (lr.abs() > 0.5).sum(),
            }
        )
agreement = pd.DataFrame(rows)
agreement.to_csv(f"{out_dir}/agreement.tsv", sep="\t", index=False)
print(f"\n## Agreement with {REF} (subfamilies with >= 50 counts)")
print(agreement.to_string(index=False))

# Young L1/IAP: counts and library-size-normalised log2 ratio to the reference
young = [y for y in YOUNG if y in ref.index]
rows = []
for m in modes:
    for s in samples:
        for y in young:
            rows.append(
                {
                    "mode": m,
                    "sample": s,
                    "subfamily": y,
                    "count": multi[m].loc[y, s],
                    "unique_count": unique[m].loc[y, s],
                    "cpm": cpm(multi[m])[s][y],
                    f"log2_vs_{REF}": np.log2(cpm(multi[m])[s][y] / cpm(ref)[s][y]),
                }
            )
young_df = pd.DataFrame(rows)
young_df.to_csv(f"{out_dir}/young_subfamilies.tsv", sep="\t", index=False)
print(f"\n## Young L1/IAP: log2(CPM / {REF} CPM)")
print(
    young_df.pivot_table(index=["subfamily"], columns=["mode"], values=f"log2_vs_{REF}")
    .loc[young, [m for m in modes if m != REF]]
    .round(2)
    .to_string()
)
print(f"\n## Young L1/IAP: uniquely counted fraction ({REF})")
u = young_df[young_df["mode"] == REF]
print(
    (
        u.set_index(["subfamily", "sample"])["unique_count"]
        / u.set_index(["subfamily", "sample"])["count"]
    )
    .unstack()
    .loc[young]
    .round(2)
    .to_string()
)

# Subfamilies that change most against the reference
rows = []
for m in modes:
    if m == REF:
        continue
    keep = ref.min(axis=1) >= 200
    lr = np.log2(cpm(multi[m]).loc[keep] / cpm(ref).loc[keep]).mean(axis=1)
    for name, v in lr.reindex(lr.abs().sort_values(ascending=False).index)[:10].items():
        rows.append({"mode": m, "subfamily": name, "mean_log2_vs_ref": round(v, 2)})
print(f"\n## Largest changes vs {REF} (>= 200 counts in all samples)")
print(pd.DataFrame(rows).to_string(index=False))
