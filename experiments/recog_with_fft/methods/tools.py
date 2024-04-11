"""
Some essential and useful functions for the algorithm behind prot-fin
"""

from typing import List, Dict, Tuple, Generator, TextIO, _GenericAlias
from sys import stderr
from tqdm import tqdm
import numpy as np
import re

# type aliases
Hash = int
WindowIndex = np.float64
ProteinID = str
Score = int
JSI = float
HashOccurence = Tuple[WindowIndex, ProteinID]
Hashes = Dict[Hash, HashOccurence]
Scores = List[Tuple[ProteinID, Tuple[WindowIndex, Score, JSI]]]
Database = Dict[Hash, List[HashOccurence]]
ProteinLookup = Dict[ProteinID, Tuple[str, int]]
ConstellationMap = List[Tuple[WindowIndex, float, float]]


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
    def __init__(self, file_name: str):
        """
        Parameters
        ----------
        file_name : str
            The name of the FASTA formatted file
        """
        with open(file_name) as f:
            self.file_name = file_name
            self.protein_count = count_appearances_in_file("^>", f)

            # validate ... TODO

    def __len__(self):
        return self.protein_count

    def __iter__(self) -> Generator[Tuple[ProteinID, str, str], None, None]:
        with open(self.file_name) as f:

            # create a progress bar and iterate over the FASTA file
            for _ in tqdm(range(len(self))):

                # find line of sequence description
                while True:
                    prot_desc = f.readline()
                    if prot_desc[0] == ">":
                        break

                # extract information from describing line
                prot_id, _, description = prot_desc.split(" ", 2)
                seq = f.readline()
                if seq[-1] == "\n":
                    seq = seq[:-1]

                # yield the extracted values, remove '>' from identifier and
                # '\n' from description
                yield prot_id[1:], description[:-1], seq


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
