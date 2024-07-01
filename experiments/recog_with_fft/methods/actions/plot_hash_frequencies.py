from tools import *
from .algorithm import get_aa_vector, create_constellation, create_hashes
import numpy as np
from matplotlib import pyplot as plt
from multiprocessing import Pool


def plot_hash_frequencies(fasta: str, outfile: str, cpu_count=1):
    fasta = Fasta(fasta)

    if cpu_count > 1:
        with Pool(cpu_count - 1) as p:
            subprocesses = p.map_async(_process, ((fasta, slice(i, None, cpu_count)) for i in range(1, cpu_count)))
            position_counts = _process((fasta, slice(0, None, cpu_count)))

            for sub_counts in subprocesses.get():
                for k, v in sub_counts.items():
                    position_counts[k] = position_counts.get(k, 0) + v

    else:
        position_counts = _process((fasta, slice(None)))

    print("Hash_Frequency,Count")
    for item in position_counts.items():
        print(*item, sep=",")

    plt.bar(position_counts.keys(), np.array(list(position_counts.values())) / len(fasta))
    plt.xlabel("Frequency of Hash in Sequence")
    plt.ylabel("Hashes")
    plt.title("Counts of average calculated Hashes")
    plt.savefig(outfile, bbox_inches='tight')


def _process(args):
    fasta, slc = args
    position_counts = {}
    for prot_id, _, seq in fasta[slc]:
        for kf in range(10):
            hashes, pos_counts = create_hashes(create_constellation(get_aa_vector(seq, kf)), prot_id, kf)
            for pos in pos_counts.values():
                position_counts[pos] = position_counts.get(pos, 0) + 1

    return position_counts
