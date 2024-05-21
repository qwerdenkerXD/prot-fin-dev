from scipy.signal import stft
from tools import *
from actions.algorithm.constellation import WINDOW_SIZE, WINDOW_TYPE
from actions.algorithm import get_aa_vector
import argparse
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
import pickle
from numpy import random, log10

random.seed(WINDOW_SIZE)


def main():
    parser = argparse.ArgumentParser(
        prog='sample-exp',
        description=__doc__,
    )
    parser.set_defaults(func=lambda _: parser.print_help())
    sub_commands = parser.add_subparsers()

    # uniref_sample_exp.py run [-s] <fasta-file> <out-file>
    run_exp_parser = sub_commands.add_parser("run", help="Run experiment")
    run_exp_parser.add_argument("fasta-file", help="The file with the protein sequences to be sampled")
    run_exp_parser.add_argument("out-file", help="The file the amplitudes will be pickled to")
    run_exp_parser.add_argument("-s", help="Number of selected samples per sequence", type=int, default=1)
    run_exp_parser.set_defaults(func=lambda args: run(getattr(args, "fasta-file"), getattr(args, "out-file"), args.s))

    # uniref_sample_exp.py plot <amplitudes-pickle>
    run_exp_parser = sub_commands.add_parser("plot", help="Plot amplitudes")
    run_exp_parser.add_argument("amplitudes-pickle", help="The file the amplitudes were pickled to")
    run_exp_parser.add_argument("out-file", help="png file for the plot output")
    run_exp_parser.set_defaults(func=lambda args: plot(getattr(args, "amplitudes-pickle"), getattr(args, "out-file")))

    args = parser.parse_args()
    args.func(args)


def run(fasta: str, out_file: str, samples_per_prot: int):
    ampl_counts = tuple({} for _ in range(WINDOW_SIZE // 2 + 1))

    for i, sample in enumerate(samples(fasta, samples_per_prot)):
        if not i % 1_000_000:
            with open(out_file, "wb") as f:
                pickle.dump(ampl_counts, f, pickle.HIGHEST_PROTOCOL)

        _, _, Zxx = stft(get_aa_vector(sample, ignore_warnings=True), nperseg=WINDOW_SIZE, noverlap=0, window=WINDOW_TYPE)

        spectrum = abs(Zxx.T)[0]
        for i, ampl in enumerate(spectrum):
            ampl = round(float(ampl), 5)
            ampl_counts[i][ampl] = ampl_counts[i].get(ampl, 0) + 1

    with open(out_file, "wb") as f:
        pickle.dump(ampl_counts, f, pickle.HIGHEST_PROTOCOL)


def plot(db: str, out_file: str):
    import re
    winsize = re.findall("WINSIZE_(\d+)", db)

    with open(db, "rb") as f:
        ampl_counts = pickle.load(f)

    plt.rcParams["font.size"] = len(ampl_counts) * 5
    fig, ax = plt.subplots()
    fig.set_figheight(len(ampl_counts) * 5)
    fig.set_figwidth(len(ampl_counts) * 5)

    raincloud = WeightedRainCloudPlot(
        range(len(ampl_counts)),
        [i.keys() for i in ampl_counts],
        [i.values() for i in ampl_counts],
        xlabel="Amplitudes",
        ylabel="Frequencies",
        winsize=winsize[0] if winsize else "")

    raincloud.plot(ax)
    plt.title("Selected Amplitudes in UniRef90", weight="bold")
    plt.savefig(out_file, bbox_inches='tight')


def samples(fasta: str, samples_per_prot: int):
    for _, _, seq in Fasta(fasta, check=False):
        sample_range = len(seq) - WINDOW_SIZE + 1
        if sample_range > 0:
            n_samples = sample_range if sample_range < samples_per_prot else samples_per_prot

            for win_start in random.choice(range(sample_range), n_samples, replace=False):
                yield seq[win_start:win_start + WINDOW_SIZE]


class WeightedRainCloudPlot:
    def __init__(self, groups: List[str], group_values: List[List[float]], group_val_weights: List[List[int]], **kwargs):
        self.groups = groups
        self.group_values = group_values
        self.group_val_weights = group_val_weights
        self.kwargs = kwargs
        self.global_max_weight = max(max(i) for i in group_val_weights)
        self.global_max_value = max(max(i) for i in group_values)

    def plot(self, ax):
        box_width = .1
        cloud_color = "gray"
        drop_color = "blue"
        max_amplitudes = []

        def get_pos(y):
            return (box_width + 1.5) * y

        print("Window_Size,Sample_Count_Mio,Frequency,5%%-Quantile,95%%-Quantile,First_Lower_Outlier,First_Upper_Outlier")

        for i, (values, val_weights) in enumerate(zip(self.group_values, self.group_val_weights)):
            pos = get_pos(i)
            box_pos = pos - box_width / 2
            max_value = max(values)
            max_amplitudes.append(max_value)
            sample_count = sum(val_weights)
            alpha = .05

            values = np.array(list(values))
            sort_by_values = np.argsort(values)
            values = values[sort_by_values]
            val_weights = np.array(list(val_weights))[sort_by_values]

            f = np.cumsum(val_weights)
            n05, first_quart, median, third_quart, n95 = \
                np.round(values[np.searchsorted(f, np.array([alpha, .25, .5, .75, 1 - alpha]) * f[-1])], 5)

            below_q1 = values[values < first_quart]
            above_q3 = values[values > third_quart]

            val_weights = val_weights / self.global_max_weight

            interquartile_range = third_quart - first_quart
            lower_whisker = first_quart - 1.5 * interquartile_range
            upper_whisker = third_quart + 1.5 * interquartile_range

            between_whiskers = (values <= upper_whisker) & (values >= lower_whisker)

            density, bin_edges = np.histogram(values[between_whiskers], weights=val_weights[between_whiskers], bins=1_000, density=True)
            raindrops = val_weights / -2 + pos - box_width
            ax.scatter(values, raindrops, s=1, color=drop_color)

            ax.fill_between(bin_edges[:-1], density / density.max() + pos, pos)

            lower_fliers = below_q1[below_q1 < lower_whisker]
            if not len(lower_fliers):
                if len(below_q1):
                    lower_whisker = below_q1.min()
                else:
                    lower_whisker = first_quart
            else:
                ax.fill_between((0, n05), pos + 1, pos - box_width - .5, color="green", alpha=.1)
                ax.text(0, pos + 1, "5%% quantile: %s" % n05, verticalalignment="top")

            upper_fliers = above_q3[above_q3 > upper_whisker]
            if not len(upper_fliers):
                if len(above_q3):
                    upper_whisker = above_q3.max()
                else:
                    upper_whisker = third_quart
            else:
                ax.fill_between((n95, self.global_max_value), pos + 1, pos - box_width - .5, color="green", alpha=.1)
                ax.text(self.global_max_value, pos + 1, "95%% quantile: %s" % n95, horizontalalignment="right", verticalalignment="top", fontsize=plt.rcParams["font.size"] / 2)

            ax.bxp(
                [{
                    'med': None,  # manually drawn median
                    'q1': first_quart,
                    'q3': third_quart,
                    'fliers': np.concatenate((lower_fliers, upper_fliers)),
                    "whislo": lower_whisker,
                    "whishi": upper_whisker
                }],
                widths=[box_width],
                positions=[box_pos],
                vert=False,
                medianprops={"linewidth": 2},
                boxprops={"linewidth": 0}  # manually drawn box
            )

            print(
                self.kwargs.get("winsize", ""),
                sample_count // 1_000_000,
                i,
                n05,
                n95,
                round(lower_fliers[-1], 5) if len(lower_fliers) else "",
                round(upper_fliers[0], 5) if len(upper_fliers) else "",
                sep=",")

            # draw box
            ax.add_patch(Rectangle((first_quart, box_pos - box_width / 4), height=box_width / 2, width=interquartile_range, color="black"))

            # draw median
            ax.add_patch(Ellipse((median, box_pos), height=box_width, width=box_width/(len(self.groups)), color="black"))

        yticks = [get_pos(i) for i in range(len(self.groups))]
        ax.set_yticks(yticks, labels=self.groups)
        max_amplitudes_axis = ax.twinx()
        max_amplitudes_axis.set_ylim(ax.get_ylim())
        max_amplitudes_axis.set_yticks(yticks, labels=max_amplitudes)
        max_amplitudes_axis.set_ylabel("Max. Amplitudes", rotation=-90, weight="bold")

        ax.set_ylabel(self.kwargs.get("ylabel"), weight="bold")
        ax.set_xlabel(self.kwargs.get("xlabel"), weight="bold")


if __name__ == '__main__':
    main()
