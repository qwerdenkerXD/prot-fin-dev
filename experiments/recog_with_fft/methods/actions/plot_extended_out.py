from tools import *
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd


def plot_extended_out(ext_out: str, plot_out: str):
    plt.figure(figsize=(50,5))

    for i, sample in enumerate(pd_read_chunkwise(ext_out), 1):
        if sample.size:
            score = sample[["JSI", "Score"]].apply(lambda x: x[0] * x[1], axis=1)
            score /= score.max()
            plt.boxplot(score, positions=[i], showfliers=False, widths=.8)

            input_fams = tuple(map(lambda x: x.split(".", 1)[0], sample["Input_Family"].iloc[0].split("|")))

            def same_fam(other):
                other_fams = tuple(map(lambda x: x.split(".", 1)[0], other.split("|")))
                return any(input_fam in other_fams for input_fam in input_fams)

            same_family = sample["Match_Family"].apply(same_fam)
            same = score[same_family]
            not_same = score[~same_family]
            plt.scatter([i] * len(same), same, marker=1, color="green")
            plt.scatter([i] * len(not_same), not_same, marker=0, color="red")

    xmin, xmax = plt.xlim()
    plt.xticks(range(int(xmin), int(xmax), 5), range(int(xmin), int(xmax), 5))
    plt.xlim(xmin - .5, xmax + .5)
    plt.xlabel("Sample index")
    plt.ylabel("JSI $\\cdot$ Score")
    plt.legend(loc="upper left", handles=(Patch(color="red", label="match not in family"), Patch(color="green", label="match in family")))
    plt.title("Found matches for different samples")
    plt.savefig(plot_out, bbox_inches='tight')
