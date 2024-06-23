"""
Some essential and useful functions for the algorithm behind prot-fin
"""

from typing import List, Dict, Tuple, Generator, TextIO, _GenericAlias
import pandas as pd
from sys import stderr
from tqdm import tqdm
import numpy as np
import re

# type aliases
Hash = int
WindowIndex = int
ProteinID = str
Score = int
JSI = float
HashOccurence = Tuple[WindowIndex, ProteinID]
Hashes = Dict[Hash, HashOccurence]
Scores = List[Tuple[ProteinID, Tuple[WindowIndex, Score, JSI]]]
Database = Dict[Hash, List[HashOccurence]]
ProteinLookup = Dict[ProteinID, Tuple[int, int]]
ConstellationMap = List[Tuple[Tuple[int, float], ...]]


def pd_read_chunkwise(csv_file: str, chunksize=10_000) -> Generator[pd.DataFrame, None, None]:
    data = pd.DataFrame()

    with open(csv_file, "r") as f:
        sample_count = count_appearances_in_file("^,,", f)

    proc_bar = tqdm(total=sample_count)

    for chunk in pd.read_csv(csv_file, sep=",", chunksize=chunksize):
        data = pd.concat([data, chunk], ignore_index=True)
        sample_ends = np.concatenate([[0], np.where(data["Rank"].isnull())], axis=None)

        if len(sample_ends) > 1:
            for i, end in enumerate(sample_ends[1:]):
                proc_bar.update(1)
                # print(i, end, data)
                sample = data.loc[sample_ends[i]: end - 1, :]
                data = data.loc[end + 1:, :]
                yield sample


def verify_type(var, ty) -> bool:
    if isinstance(ty, type):
        assert ty not in (dict, list, tuple), \
            "verify_type: Use types from typing instead of builtins for dict, list, tuple - got: " + ty.__name__
        return type(var) is ty

    if ty.__origin__ is list:
        if not isinstance(var, list):
            return False
        if False in [verify_type(item, ty.__args__[0]) for item in var]:
            return False

    elif ty.__origin__ is dict:
        if not isinstance(var, dict):
            return False
        if False in [verify_type(item, ty.__args__[0]) for item in var.keys()]:
            return False
        if False in [verify_type(item, ty.__args__[1]) for item in var.values()]:
            return False

    elif ty.__origin__ is tuple:
        if not isinstance(var, tuple):
            return False
        if ty.__args__[-1] is ...:
            if False in [verify_type(item, ty.__args__[0]) for item in var]:
                return False

        else:
            if len(var) != len(ty.__args__):
                return False
            if False in [verify_type(item, type_) for item, type_ in zip(var, ty.__args__)]:
                return False

    return True


class Fasta:
    """
    A class used for convenient iteration over a FASTA file's contents.

    ...

    Attributes
    ----------
    file_name : str
        The name of the FASTA formatted file
    protein_count : int
        The number of sequences stored in the FASTA file
    """
    def __init__(self, file_name: str, check=True):
        """
        Parameters
        ----------
        file_name : str
            The name of the FASTA formatted file
        """
        with open(file_name) as f:
            self.file_name = file_name
            if check:
                self.protein_count = count_appearances_in_file("^>", f)

                # validate ... TODO or validate during iteration like (re.match("^[A-Z]+$", seq) is not None)

            else:
                self.protein_count = None

    def __len__(self):
        if self.protein_count is None:
            raise TypeError("Can't predict protein count of unchecked Fasta")
        return self.protein_count

    def __iter__(self) -> Generator[Tuple[ProteinID, str, str], None, None]:
        return self[:]

    def __getitem__(self, key) -> Generator[Tuple[ProteinID, str, str], None, None]:
        assert type(key) is slice
        if self.protein_count is None:
            start, _, step = key.indices(0)
            stop = key.stop

            def open_range(start, step=1):
                while True:
                    yield start
                    start += step

            pbar = tqdm(open_range(start, step))
        else:
            start, stop, step = key.indices(self.protein_count)
            pbar = tqdm(range(start, stop, step))
        assert step > 0, "negative steps currently not supported"

        with open(self.file_name) as f:
            # create a progress bar and iterate over the FASTA file
            current_prot = 0
            for i in pbar:

                # find line of sequence description
                while current_prot <= i and (prot_desc := f.readline()):
                    current_prot += prot_desc[0] == ">"
                if not prot_desc:
                    break  # because EOF

                # read sequence
                seq = ""
                while (seq_segment := f.readline()) and seq_segment[0] != ">":
                    pointer = f.tell()
                    seq += seq_segment
                seq = seq.replace("\n", "")

                # go to last line -> readline returns header of next protein
                f.seek(pointer)

                # extract information from describing line
                header = prot_desc.split(" ", 1)
                if len(header) == 2:
                    prot_id, description = header
                else:
                    prot_id, description = header[0][:-1], "\n"

                # yield the extracted values, remove '>' from identifier and
                # '\n' from description
                yield prot_id[1:], description[:-1], seq

        if stop is not None:
            assert current_prot == max(start, stop - (stop - start - 1) % step)


def count_appearances_in_file(pattern, file: TextIO):
    count = 0
    file.seek(0)
    while (buffer := file.read(1024 ** 2)):
        count += len(re.findall(pattern, buffer, re.MULTILINE))
    file.seek(0)

    return count


def eprint(*args, **kwargs):
    print(*args, file=stderr, **kwargs)


def warn(*args, **kwargs):
    eprint("WARNING:", *args, **kwargs)
