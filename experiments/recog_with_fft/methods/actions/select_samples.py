from typing import List, Dict
from tools import Fasta, ProteinID
import pandas as pd
from numpy import random

Samples = Dict[ProteinID, str]

random.seed(1)


def select_samples(mapman: str, protein_file: str, samples_per_family: int):
    """
    Analyzes the mapman bins to select random functional different
    proteins as train data.

    The data is printed to stdout

    ...

    Parameters
    ----------
    mapman : str
        Path to the mapman result for the proteins.
    protein_file : str
        Path to the FASTA formatted file storing the amino acid sequences for
        the proteins used for mapman.
    samples_per_family : int
        Number of samples that will be selected per each family.
    """

    # read the mapman reference
    map_data: pd.Series = read_mapman(mapman)
    families = map_data.value_counts()
    samples: Samples = get_samples(families, samples_per_family, map_data)

    print_samples(samples, protein_file)


def read_mapman(mapman: str) -> pd.Series:
    data = pd.read_csv(mapman, sep="\t", quotechar="'", index_col=1, usecols=["IDENTIFIER", "BINCODE"]).squeeze()
    data = data[data.index.notna() & ~data.str.startswith("50.") & ~data.str.startswith("35.")]
    return data


def print_samples(samples, protein_file):
    for prot_id, _, seq in Fasta(protein_file):
        if prot_id.lower() in samples:
            print(">%s  %s\n%s" % (prot_id, samples[prot_id.lower()], seq))


def get_samples(families: pd.Series, samples_per_family: int, map_data: pd.DataFrame) -> Samples:
    selected_samples: Samples = {}
    families = families[families > 1].groupby(lambda bincode: int(bincode.split(".", 1)[0]))

    for _, fam in families:
        chosen_idx = random.choice(len(fam), samples_per_family, replace=False)
        assert len(fam.index[chosen_idx]) == samples_per_family
        for fam_id in fam.index[chosen_idx]:
            members = map_data[(map_data == fam_id)]

            members = members[~members.index.isin(selected_samples.keys())]
            prot_id = members.index[random.choice(len(members))]
            selected_samples[prot_id] = fam_id

    return selected_samples
