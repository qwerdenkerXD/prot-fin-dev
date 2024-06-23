from tools import *

FREQUENCY_BITS = 5
DIFFERENCE_BITS = 3


def create_hashes(
        constellation_map: ConstellationMap,
        prot_id: ProteinID,
        kidera_factor: int
        ) -> Hashes:
    """
    Creates combinatorial Hashes from a constellation map for efficient
    database searches

    ...

    Parameters
    ----------
    constellation_map : ConstellationMap
        A list of coordinates, index-frequency pairs, the result of the STFT
        in create_constellation. It is assumed as pre-sorted by index
    prot_id : ProteinID
        The identifier for the protein the hashes are generated for. If it is
        unknown, just pass a custom

    Returns
    -------
    A dictionary of hashes pointing to the index of their occurence
    """

    hashes: Hashes = {}

    # Iterate through the constellation map
    for idx, freqs in enumerate(constellation_map):
        occ = (idx, prot_id)
        # Iterate through the next pairs to produce combinatorial hashes
        for freq, _, quantile in freqs:
            hash_count = len(hashes)
            for diff, other_freqs in enumerate(constellation_map[idx + 1:idx + 2**DIFFERENCE_BITS]):
                for other_freq, _, other_quantile in other_freqs:
                    # Produce a 32 bit hash
                    hash_: Hash = create_hash((kidera_factor, 4), (quantile, 1), (other_quantile, 1), (diff, DIFFERENCE_BITS), (other_freq, FREQUENCY_BITS), (freq, FREQUENCY_BITS))
                    hashes[hash_] = occ
            if len(hashes) == hash_count:  # -> no hashes created for the quantile's frequency -> combining with foo frequency
                other_freq = 2 ** FREQUENCY_BITS - 1
                other_quantile = 0
                diff = 0
                hash_: Hash = create_hash((kidera_factor, 4), (quantile, 1), (other_quantile, 1), (diff, DIFFERENCE_BITS), (other_freq, FREQUENCY_BITS), (freq, FREQUENCY_BITS))
                hashes[hash_] = occ
    return hashes


def create_hash(*args) -> Hash:
    hash_: Hash = 0
    bits = 0
    for val, shift in args:
        assert int(val) < 2 ** shift, "%s too big for %s bit" % (val, shift)

        # move bits to the left to make space for the next value
        hash_ <<= shift

        # insert the value's bits
        hash_ |= int(val)

        bits += shift
    assert bits <= 32, "Hash exceeds 32 bit"

    return hash_
