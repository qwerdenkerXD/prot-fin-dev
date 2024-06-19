from typing import Dict, List, Tuple, TextIO
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

    mapman_data = pd.read_csv(mapman, sep="\t", quotechar="'", index_col=1, usecols=["IDENTIFIER", "BINCODE"])
    mapman_data = mapman_data[mapman_data.index.notna()]
    mapman_data = mapman_data["BINCODE"].apply(str)

    # data will be collected into a dataframe
    evaluation = pd.DataFrame(columns=["Sample_ID", "First_Match_Count", "Sample_In_First_Matches", "Sequence_Length", "Sample_Hashes", "Sample_JSI", "Sample_Score", "Top_JSI", "Top_Score", "F1_Score"])
    for matches in pd_read_chunkwise(protfin_out_file):
        if matches.size:
            input_sample = matches["Input_Protein_ID"].iloc[0]
            sample_result = matches[matches["Match_Protein_ID"] == input_sample]
            # Q94A34 Q9ZSE4
            input_fams = np.array(mapman_data.loc[input_sample.lower()], ndmin=1)

            def same_fam(other):
                return other in input_fams

            positives = matches["Rank"] <= int(len(matches.index) * .05) + 1
            pos_matches_fams = mapman_data[matches["Match_Protein_ID"][positives].apply(str.lower)]
            true_positives = pos_matches_fams.apply(same_fam).groupby(pos_matches_fams.index).max().sum()
            precision = true_positives / positives.sum()

            neg_matches_fams = mapman_data[matches["Match_Protein_ID"][~positives].apply(str.lower)]
            false_negatives = neg_matches_fams.apply(same_fam).groupby(neg_matches_fams.index).max().sum()
            recall = true_positives / (true_positives + false_negatives) if true_positives + false_negatives else 0
            f1_score = (2 * precision * recall) / (precision + recall) if precision + recall else 0

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
                f1_score
            )

    # write to stdout
    print(evaluation.to_csv(index=False, float_format="%g"), end="")
