from tools import *
import pandas as pd
import numpy as np
from os import environ as env

KIDERA_FACTOR = ["Helix/bend preference", "Side-chain size", "Extended structure preference", "Hydrophobicity", "Double-bend preference", "Partial specific volume", "Flat extended preference", "Occurrence in alpha region", "pK-C", "Surrounding hydrophobicity"]\
    .index("Hydrophobicity")

KIDERA_TABLE = pd.DataFrame.from_dict({
    "A": [-1.56, -1.67, -0.97, -0.27, -0.93, -0.78, -0.2,  -0.08,  0.21, -0.48 ],
    "C": [ 0.12, -0.89,  0.45, -1.05, -0.71,  2.41,  1.52, -0.69,  1.13,  1.1  ],
    "D": [ 0.58, -0.22, -1.58,  0.81, -0.92,  0.15, -1.52,  0.47,  0.76,  0.7  ],
    "E": [-1.45,  0.19, -1.61,  1.17, -1.31,  0.4,   0.04,  0.38, -0.35, -0.12 ],
    "F": [-0.21,  0.98, -0.36, -1.43,  0.22, -0.81,  0.67,  1.1,   1.71, -0.44 ],
    "G": [ 1.46, -1.96, -0.23, -0.16,  0.1,  -0.11,  1.32,  2.36, -1.66,  0.46 ],
    "H": [-0.41,  0.52, -0.28,  0.28,  1.61,  1.01, -1.85,  0.47,  1.13,  1.63 ],
    "I": [-0.73, -0.16,  1.79, -0.77, -0.54,  0.03, -0.83,  0.51,  0.66, -1.78 ],
    "K": [-0.34,  0.82, -0.23,  1.7,   1.54, -1.62,  1.15, -0.08, -0.48,  0.6  ],
    "L": [-1.04,  0.0,  -0.24, -1.1,  -0.55, -2.05,  0.96, -0.76,  0.45,  0.93 ],
    "M": [-1.4,   0.18, -0.42, -0.73,  2.0,   1.52,  0.26,  0.11, -1.27,  0.27 ],
    "N": [ 1.14, -0.07, -0.12,  0.81,  0.18,  0.37, -0.09,  1.23,  1.1,  -1.73 ],
    "P": [ 2.06, -0.33, -1.15, -0.75,  0.88, -0.45,  0.3,  -2.3,   0.74, -0.28 ],
    "Q": [-0.47,  0.24,  0.07,  1.1,   1.1,   0.59,  0.84, -0.71, -0.03, -2.33 ],
    "R": [ 0.22,  1.27,  1.37,  1.87, -1.7,   0.46,  0.92, -0.39,  0.23,  0.93 ],
    "S": [ 0.81, -1.08,  0.16,  0.42, -0.21, -0.43, -1.89, -1.15, -0.97, -0.23 ],
    "T": [ 0.26, -0.7,   1.21,  0.63, -0.1,   0.21,  0.24, -1.15, -0.56,  0.19 ],
    "V": [-0.74, -0.71,  2.04, -0.4,   0.5,  -0.81, -1.07,  0.06, -0.46,  0.65 ],
    "W": [ 0.3,   2.1,  -0.72, -1.57, -1.16,  0.57, -0.48, -0.4,  -2.3,  -0.6  ],
    "Y": [ 1.38,  1.48,  0.8,  -0.56,  0.0,  -0.68, -0.31,  1.03, -0.05,  0.53 ]
})

KIDERA_MIN = abs(KIDERA_TABLE.to_numpy().min())


def get_aa_vector(
        seq: str,
        factor=KIDERA_FACTOR,
        normalize=True,
        ignore_warnings=False
        ) -> np.ndarray:
    """
    Transform an amino acid sequence into a vector of floats from the
    Kidera factor table for the selected factor

    ...

    Parameters
    ----------
    seq : str
        The amino acid sequence to be transformed
    factor : int
        The index of the Kidera factor

    Returns
    -------
    A numpy array of 32-bit floats
    """

    # select the specified Kidera factor from the table
    sel_factor: pd.Series = KIDERA_TABLE.iloc[factor]

    # normalizing to non-negatives by adding the absolute of the global minimum
    if normalize:
        sel_factor = KIDERA_TABLE.iloc[factor] + KIDERA_MIN + 1

    # define symbols representing multiple amino acids
    extend_selected_factor(sel_factor)

    # transform the sequence
    transformed_seq = transform_seq(seq, sel_factor, ignore_warnings)

    return np.array(transformed_seq)


def extend_selected_factor(sel_factor: pd.Series):
    special_aa = {
        "X": sel_factor.keys(),  # any aminoacid
        "B": ["D", "N"],
        "Z": ["E", "Q"],
        "J": ["I", "L"],
        "Ψ": ["I", "L", "M", "V"],
        "Ω": ["F", "W", "Y", "H"],
        "Φ": ["I", "L", "M", "V", "F", "W", "Y"],
        "ζ": ["D", "E", "H", "K", "N", "Q", "R", "S", "T"],
        "Π": ["A", "G", "P", "S"],
        "+": ["K", "R", "H"],
        "-": ["D", "E"]
    }

    # extend the factor data with the multi-representing symbols
    for s in special_aa:
        sel_factor[s] = sel_factor[special_aa[s]].mean()


def transform_seq(seq, sel_factor, ignore_warnings=False) -> List[np.float32]:
    transformed_seq: List[np.float32] = []
    for aa in seq:
        value = sel_factor.get(aa)

        # currently 'O' and 'U' are unknown
        if value is None:
            value = 0
            if not ignore_warnings:
                warn(f"No known values for Kidera factors for {aa} -> treating as zero")
        transformed_seq.append(np.float32(value))

    return transformed_seq
