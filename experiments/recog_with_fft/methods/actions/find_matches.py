from typing import Dict, List, Tuple
import pandas as pd
from tools import *
from .algorithm import hashes_from_seq
import pickle
import sys
import os

Matches = List[Tuple[WindowIndex, WindowIndex]]
ScoresByOffset = Dict[WindowIndex, Score]
MatchesPerProt = Dict[ProteinID, ScoresByOffset]
ScoresMap = Dict[ProteinID, Tuple[WindowIndex, Score, JSI]]

COLUMNS = [
    "Rank",
    "Match_Protein_ID",
    "JSI",
    "Score",
    "Input_Protein_ID",
    "Input_Sequence_Length",
    "Input_Found_Hashes"
]


def find_matches(
        fasta_file: str,
        db_in: str,
        filter_quantile=1.0
        ):
    """
    Find matches for the proteins defined in the FASTA file
    and print them to stdout

    ...

    Parameters
    ----------
    fasta_file : str
        The path to the FASTA formatted file containing all protein sequences
        of interest
    db_in : str
        Name of the file storing the trained database
    filter_quantile : float
        Quantile of hashes to be kept in database
    """

    assert filter_quantile > 0 and filter_quantile <= 1

    (database, protein_lookup), hash_blacklist = get_filtered_db(db_in, filter_quantile)

    print(*COLUMNS, sep=",")

    # print the matches with description and score
    for input_id, _, seq in Fasta(fasta_file):
        # create the combinatorial hashes for the sequence
        hashes: Hashes = hashes_from_seq(seq, None)
        hashes = {k: v for k, v in hashes.items() if k not in hash_blacklist}

        # calculate the scores for proteins in the database
        scored_matches: ScoresMap = score_prots(hashes, database, protein_lookup)

        result = get_result_frame(scored_matches, input_id, len(seq), len(hashes))

        result.loc[len(result.index)] = None
        print(result.sort_values("Rank").to_csv(index=False, header=False, float_format="%g"), end="")


def get_filtered_db(db_in: str, filter_quantile: float) -> Tuple[Tuple[Database, ProteinLookup], List[Hash]]:
    # load databases
    with open(db_in, 'rb') as f:
        database, protein_lookup = pickle.load(f)

    prev_size = os.path.getsize(db_in)

    db, hash_blacklist = filter_db(database, protein_lookup, filter_quantile)

    size_now = sys.getsizeof(pickle.dumps(db, pickle.HIGHEST_PROTOCOL))
    eprint(int(size_now / prev_size * 100), r"%", " (%.2fMB) of database size used by quantile filter" % (size_now / (1024**2)), sep="")

    return db, hash_blacklist


def filter_db(database: Database, protein_lookup: ProteinLookup, filter_quantile: float) -> Tuple[Tuple[Database, ProteinLookup], List[Hash]]:
    # filter db
    hash_frequencies = np.array(sorted(len(v) for v in database.values()))
    f = np.cumsum(hash_frequencies)
    quantile_value = hash_frequencies[np.searchsorted(f, filter_quantile * f[-1])]

    hash_blacklist = [hash_ for hash_, v in database.items() if len(v) > quantile_value]
    database = {k: v for k, v in database.items() if len(v) <= quantile_value}

    # update lookup
    hash_counts = {k: 0 for k in protein_lookup}
    for prots in database.values():
        for _, p in prots:
            hash_counts[p] += 1
    for p, val in protein_lookup.items():
        protein_lookup[p] = (val[0], hash_counts[p])

    return (database, protein_lookup), hash_blacklist


def get_result_frame(
        scored_matches: ScoresMap,
        input_id: str,
        seq_len: int,
        hash_count: int
        ) -> pd.DataFrame:

    result = pd.DataFrame(columns=COLUMNS)

    if scored_matches:
        result["Match_Protein_ID"] = scored_matches.keys()
        _, scores, jsi = zip(*scored_matches.values())
        result["JSI"] = jsi
        result["Score"] = scores
        result["Input_Protein_ID"] = input_id
        result["Input_Sequence_Length"] = seq_len
        result["Input_Found_Hashes"] = hash_count

        result["Rank"] = result[["JSI", "Score"]].apply(lambda x: x[0] * x[1], axis=1).rank(method="dense", ascending=False)

    return result


def score_prots(
        hashes: Hashes,
        database: Database,
        protein_lookup: ProteinLookup
        ) -> ScoresMap:
    """
    Scores the proteins of the database by their suitability to the given hashes

    ...

    Parameters
    ----------
    hashes : Hashes
        A dictionary of combinatorial hashes generated from a
        sample protein sequence pointing to the index of their occurence and
        the id of the protein
    database: Database
        A dictionary of combinatorial hashes generated from known protein
        sequences pointing to all proteins they occur in, including the position
    protein_lookup: ProteinLookup
        A dictionary of protein identifiers pointing to their description and
        number of combinatorial hashes they have, necessary to calclulate the
        Jaccard Similarity Score

    Returns
    -------
    A list of match protein identifiern and their respective scores.
        The float value is the Jaccard Similarity Index and describes the ratio
        of how many hashes of all possibles for a sequence are found,
        the second value is an additional score to describe how good the found
        hashes fit to the respective protein,
        the first value describes the offset of the matching hashes to their
        original position in the protein sequence
    """

    matches_per_prot: MatchesPerProt = get_matches_per_prot(hashes, database)

    # stores all identifiers of proteins that have matching hashes, pointing
    # to the distance to the match region in their sequence, the score and the
    # Jaccard Similarity Index
    scores_map: ScoresMap = {}

    # calculate the scores for each protein based on the matching hashes
    for match_prot_id, offsets in matches_per_prot.items():

        # calculate the jaccard similarity index as one scoring value
        sample_cardinality: int = len(hashes)
        _, match_cardinality = protein_lookup[match_prot_id]

        intersection_cardinality: int = sum(offsets.values())
        union_cardinality: int = sample_cardinality + match_cardinality - intersection_cardinality

        # union cardinality can't be zero, because:
        # sample_cardinality = 0 -> no found hashes -> no matches per prot
        # match_cardinality = 0 -> no found hashes -> not in database -> not in matches per prot
        # -> union_cardinality >= 2
        jacc_sim_index: JSI = intersection_cardinality / union_cardinality

        # find the most related offset, as it describes the biggest matching
        # constellation of found hashes for a protein
        max_offset, max_frequency = get_max_offset(offsets)

        # take the frequency of the most related offset as score
        scores_map[match_prot_id] = (max_offset, max_frequency, jacc_sim_index)

    return scores_map


def get_max_offset(prot_scores_by_offset: ScoresByOffset) -> Tuple[WindowIndex, Score]:
    max_ = (0, 0)
    for offset, frequency in prot_scores_by_offset.items():
        if frequency > max_[1]:
            max_ = (offset, frequency)

    return max_


def get_matches_per_prot(hashes: Hashes, database: Database) -> MatchesPerProt:
    # stores all found hashes that exist for a known protein
    matches_per_prot: MatchesPerProt = {}

    # find the matches per protein
    for hash_, (sample_index, _) in hashes.items():
        if hash_ in database:
            matching_occurences: List[HashOccurence] = database[hash_]

            # for each known protein for the hash, keep the indices of its
            # occurrence in the sample and the known sequence
            for source_index, match_prot_id in matching_occurences:
                if (offsets := matches_per_prot.get(match_prot_id)) is None:
                    offsets = {}
                    matches_per_prot[match_prot_id] = offsets

                offset = sample_index - source_index
                offsets[offset] = offsets.get(offset, 0) + 1

    return matches_per_prot
