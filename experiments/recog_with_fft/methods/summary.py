import pandas as pd
from glob import glob
from sys import argv as args
import re

res = pd.DataFrame(columns=[
    "Window_Size",
    "n_Peaks",
    "Overlap",
    "Median_Hits",
    "Average_Hits",
    "Self_Matches",
    "Unique_Self_Matches",
    "Mean_F1_Score"
    ]
)
for file in args[1:]:
    summary = pd.read_csv(file)
    res.loc[len(res.index)] = (
        *re.findall("WINSIZE_(\d+)_NPEAKS_(\d+)_OVERLAP_(\d+)", file)[0],
        summary["First_Match_Count"].median(),
        round(summary["First_Match_Count"].mean(), 2),
        (summary["Sample_In_First_Matches"]).sum(),
        (summary["Sample_In_First_Matches"] & (summary["First_Match_Count"] == 1)).sum(),
        summary["F1_Score"].mean()
    )

print(res.sort_values("Mean_F1_Score", ascending=False).to_csv(index=False, float_format="%g"))
