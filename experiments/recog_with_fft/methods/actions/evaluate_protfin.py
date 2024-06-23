from typing import Dict, List, Tuple, TextIO
from .select_samples import read_mapman
from tools import ProteinID, JSI, Score
from tools import pd_read_chunkwise
from tqdm import tqdm
import pandas as pd
import numpy as np
import re


Matches = Dict[ProteinID, Tuple[JSI, Score]]


def evaluate_protfin(protfin_out_file: str, mapman: str):
    """
    Summarizes the output of prot-fin's find-match output.
    For each sample protein, it collects the following information into csv:
      - sample protein id                         -> Sample_ID
      - number of found matches                   -> First_Match_Count
      - occurs the sample in the found matches    -> Sample_In_First_Matches
      - sequence length of sample                 -> Sequence_Length
      - created hashes for sample                 -> Sample_Hashes
      - Jaccard Similarity Score of sample        -> Sample_JSI
      - Score of sample                           -> Sample_Score
      - Jaccard Similarity Score of found matches -> Sample_JSI
      - Score of found matches                    -> Sample_Score

    The data is printed to stdout

    ...

    Parameters
    ----------
    protfin_out_file : str
        Path to the file the output of protfin was written to
    """

    mapman_data = read_mapman(mapman)

    # data will be collected into a dataframe
    evaluation = pd.DataFrame(columns=["Sample_ID", "First_Match_Count", "Sample_In_First_Matches", "Sequence_Length", "Sample_Hashes", "Sample_JSI", "Sample_Score", "Top_JSI", "Top_Score", "F1_Score", "Liberal_F1_Score", "Precision", "Precision_Liberal"])
    for matches in pd_read_chunkwise(protfin_out_file):
        if matches.size:
            input_sample = matches["Input_Protein_ID"].iloc[0]
            sample_result = matches[matches["Match_Protein_ID"] == input_sample]

            input_fams = np.array(mapman_data[input_sample.lower()], ndmin=1, dtype=str)

            members = mapman_data.isin(input_fams).groupby(mapman_data.index).max()
            match_prots = matches["Match_Protein_ID"].apply(str.lower)
            matches_fams = members[match_prots[match_prots.isin(members.index)]][::-1]
            member_counts = matches_fams.cumsum()

            true_positives = member_counts.iloc[-1]
            false_negatives = members.sum() - true_positives
            recall = true_positives / (true_positives + false_negatives)

            def f1_score(positives: int):
                precision = true_positives / positives
                return (2 * precision * recall) / (precision + recall), precision

            f1, f1_prec = f1_score(len(matches.index))  # normal f1_score has all found matches as positives
            lib_f1, lib_prec = f1_score(len(matches.index) - member_counts.value_counts().get(0, 0))  # liberal f1_score has all matches until last family member as positives, all negatives have cumsum value 0
            # insert the data into the dataframe
            evaluation.loc[len(evaluation.index)] = (
                input_sample,  # Sample_ID
                (matches["Rank"] == 1).sum(),  # First_Match_Count
                input_sample in matches[matches["Rank"] == 1]["Match_Protein_ID"].values,  # Sample_In_First_Matches
                matches["Input_Sequence_Length"].iloc[0],  # Sequence_Length
                matches["Input_Found_Hashes"].iloc[0],  # Sample_Hashes
                sample_result["JSI"].iloc[0] if len(sample_result) else None,  # Sample_JSI
                sample_result["Score"].iloc[0] if len(sample_result) else None,  # Sample_Score
                matches[matches["Rank"] == 1]["JSI"].max(),  # Top_JSI
                matches[matches["Rank"] == 1]["Score"].max(),  # Top_Score
                f1,
                lib_f1,
                f1_prec,
                lib_prec
            )

    # write to stdout
    print(evaluation.to_csv(index=False, float_format="%g"), end="")
