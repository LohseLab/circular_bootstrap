#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
from matplotlib import pyplot as plt


def resample(tt, xx, stat, j):
    test_l = list(tt[lab_col])
    off_x_list = [test_l[(i + int(xx[j])) % len(test_l)] for i, x in enumerate(test_l)]
    randomD = list(tt[stat][np.array(off_x_list)])
    return randomD


def resample_wrapper(df, xx, stat, chroms=None, global_test=False):
    if global_test:
        return np.mean(resample(df, xx, stat, j=0))

    else:
        pp = []
        for j, chrom in enumerate(chroms):
            tt = df[df[chr_col] == chrom]
            pp.extend(resample(tt, xx, stat, j))

        return np.mean(pp)


def plot_results(hist, obs_mean, stat, hist_label):
    fig, axs = plt.subplots(figsize=(15, 10))

    plt.hist(hist, color="grey", label=str(hist_label), bins=int(len(hist) / 50))
    plt.axvline(x=obs_mean, color="red", label="Observed", linestyle="--")

    axs.set_xlabel(str(stat))
    axs.legend(loc="upper center")

    f_prefix = str(stat).replace(" ", "_")
    fig.savefig(f"{f_prefix}.circular_bootstrap.png", dpi=200)


gIMble_f = sys.argv[1]
out_prefix = gIMble_f.split('.')[0]
gimble_wins_df = pd.read_table(gIMble_f, delimiter="\t")
chr_col = gimble_wins_df.columns[0]
lab_col = gimble_wins_df.columns[1]
stats = gimble_wins_df.columns[2:]
chroms = gimble_wins_df[chr_col].unique()

# setting the random seed
seed = 675
n = 1000
rng = np.random.default_rng(seed)
offsets = np.zeros((len(chroms), n))
for i, chrom in enumerate(chroms):
    length = np.sum(gimble_wins_df[chr_col] == chrom)
    offsets[i] = rng.integers(length, size=n)
offsets = offsets.T

global_offset = [[x] for x in rng.integers(gimble_wins_df.shape[0], size=n)]

results_chrom = {}
results_global = {}

for stat in stats:
    result_series_chrom = np.array(
        [
            resample_wrapper(
                gimble_wins_df, offsets[i], stat, chroms, global_test=False
            )
            for i in range(n)
        ]
    )

    result_series_global = np.array(
        [
            resample_wrapper(
                gimble_wins_df, global_offset[i], stat, chroms=None, global_test=True
            )
            for i in range(n)
        ]
    )

    results_chrom[stat] = result_series_chrom
    results_global[stat] = result_series_global

    plot_results(
        results_chrom[stat],
        gimble_wins_df[gimble_wins_df[lab_col] == True][stat].mean(),
        stat,
        "Circular bootstrap (per chromosome)",
    )
    plt.close()

pd.DataFrame.from_dict(results_chrom).to_csv(
    f"{out_prefix}.chrom_resamples.tsv", sep="\t", index=None
)
pd.DataFrame.from_dict(results_global).to_csv(
    f"{out_prefix}.global_resamples.tsv", sep="\t", index=None
)
