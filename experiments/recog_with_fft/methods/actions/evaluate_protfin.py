from typing import Dict, List, Tuple, TextIO
from tools import ProteinID, JSI, Score
from tools import pd_read_chunkwise
from tqdm import tqdm
import pandas as pd
import re


Matches = Dict[ProteinID, Tuple[JSI, Score]]


def evaluate_protfin(protfin_out_file: str):
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

    # data will be collected into a dataframe
    evaluation = pd.DataFrame(columns=["Sample_ID", "First_Match_Count", "Sample_In_First_Matches", "Sequence_Length", "Sample_Hashes", "Sample_JSI", "Sample_Score", "Top_JSI", "Top_Score", "Median_Score_Input_Family", "Mean_Score_Input_Family", "Median_Score_Not_Input_Family", "Mean_Score_Not_Input_Family"])
    for matches in pd_read_chunkwise(protfin_out_file):
        if matches.size:
            input_sample = matches["Input_Protein_ID"].iloc[0]
            sample_result = matches[matches["Match_Protein_ID"] == input_sample]

            input_fams = tuple(map(lambda x: x.split(".", 1)[0], matches["Input_Family"].iloc[0].split("|")))

            def same_fam(other):
                other_fams = tuple(map(lambda x: x.split(".", 1)[0], other.split("|")))
                return any(input_fam in other_fams for input_fam in input_fams)

            same_family = matches["Match_Family"].apply(same_fam)
            score = matches[["JSI", "Score"]].apply(lambda x: x[0] * x[1], axis=1)
            same = score[same_family]
            not_same = score[~same_family]

            # insert the data into the dataframe
            evaluation.loc[len(evaluation.index)] = (
                input_sample,  # Sample_ID
                (matches["Rank"] == 1).sum(),  # First_Match_Count
                input_sample in matches[matches["Rank"] == 1]["Match_Protein_ID"].values,  # Sample_In_First_Matches
                matches["Input_Sequence_Length"].iloc[0],  # Sequence_Length
                matches["Input_Found_Hashes"].iloc[0],  # Sample_Hashes
                sample_result["JSI"].iloc[0],  # Sample_JSI
                sample_result["Score"].iloc[0],  # Sample_Score
                matches[matches["Rank"] == 1]["JSI"].max(),  # Top_JSI
                matches[matches["Rank"] == 1]["Score"].max(),  # Top_Score
                same.median(),  # Median_Score_Input_Family
                same.mean(),  # Mean_Score_Input_Family
                not_same.median(),  # Median_Score_Not_Input_Family
                not_same.mean()  # Mean_Score_Not_Input_Family
            )

    # write to stdout
    print(evaluation.to_csv(index=False, float_format="%g"), end="")
