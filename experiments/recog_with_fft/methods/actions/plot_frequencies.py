from tools import *
from .algorithm import create_constellation, get_aa_vector
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as Cmap, Normalize
from matplotlib.cm import ScalarMappable
from multiprocessing import Pool
from .algorithm.constellation import WINDOW_SIZE
from scipy.fft import fftfreq

FREQS = list(filter(lambda x: x >= 0, fftfreq(WINDOW_SIZE)))


def plot_frequencies(prot_file: str, out_file: str, cpu_count=1):
    fasta = Fasta(prot_file)
    protein_count = len(fasta)

    if cpu_count > 1:
        with Pool(cpu_count - 1) as p:
            subprocesses = p.map_async(_process, ((fasta, slice(i, None, cpu_count)) for i in range(1, cpu_count)))
            sel_freqs = _process((fasta, slice(0, None, cpu_count)))
            for res_sel_freqs in subprocesses.get():
                for res_freq, (res_freq_count, res_prot_count) in res_sel_freqs.items():
                    freq_count, prot_count = sel_freqs[res_freq]
                    sel_freqs[res_freq] = (freq_count + res_freq_count, prot_count + res_prot_count)
    else:
        sel_freqs = _process((fasta, slice(None)))

    if sel_freqs:
        cmap = Cmap.from_list(
            'custom_cmap',
            [(0, 'orange'), (1/3, 'red'), (1, 'black')]
        )
        norm = Normalize(vmin=0, vmax=protein_count)

        freqs, counts = zip(*sorted(filter(lambda x: x[1][0] > 0, sel_freqs.items())))
        freq_counts, prot_counts = zip(*counts)

        fig = plt.figure(figsize=(10, 5))
        plt.scatter(freqs, freq_counts, marker='o', c=list(map(lambda c: cmap(norm(c)), prot_counts)))
        for f, f_count in zip(freqs, freq_counts):
            plt.text(
                f,
                f_count - (max(freq_counts)-min(freq_counts)) / 100 * 1.5,
                round(f, 2),
                ha="center",
                va="top"
            )
        plt.colorbar(ScalarMappable(cmap=cmap, norm=norm)).set_label("Sequences", rotation=-90, va="bottom")
        plt.xlabel("Frequencies")
        plt.ylabel("Absolute count in all sequences")
        plt.title("Occurences of selected STFT frequences")
        plt.savefig(out_file, bbox_inches='tight')


def _process(args) -> Dict[float, Tuple[int, int]]:
    fasta, prot_slice = args
    selected_freqs = dict.fromkeys(FREQS, (0, 0))

    for _, _, seq in fasta[prot_slice]:
        constellation_map = create_constellation(get_aa_vector(seq))
        all_freqs = set()
        for window in constellation_map:
            for freq_idx, _, _ in window:
                freq = FREQS[freq_idx]
                freq_count, prot_count = selected_freqs[freq]
                selected_freqs[freq] = (freq_count + 1, prot_count)
                all_freqs.add(freq)

        for freq in all_freqs:
            freq_count, prot_count = selected_freqs[freq]
            selected_freqs[freq] = (freq_count, prot_count + 1)

    return selected_freqs
