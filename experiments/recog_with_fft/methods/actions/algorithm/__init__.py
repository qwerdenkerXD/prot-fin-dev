from tools import *
from .kidera import get_aa_vector
from .constellation import create_constellation
from .hash_gen import create_hashes
from os import environ as env


def hashes_from_seq(seq: str, prot_id: str) -> Hashes:
    """
    Generate the combinatorial hashes from an amino acid sequence

    ...

    Parameters
    ----------
    seq : str
        A sequence of amino acid in their one letter codes
    prot_id : ProteinID
        The identifier for the protein of the passed sequence. If unknown,
        just pass a custom
    """

    hashes = {}

    KF = env.get("KIDERA_FACTOR")
    for kf in ((int(KF),) if KF is not None else range(10)):
        # start the pipeline
        aa_vec = get_aa_vector(seq, kf)
        constellation = create_constellation(aa_vec)
        hashes = {**hashes, **create_hashes(constellation, prot_id, kf)[0]}

    return hashes
