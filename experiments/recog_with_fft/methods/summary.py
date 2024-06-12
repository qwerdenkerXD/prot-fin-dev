import pandas as pd
from glob import glob
from sys import argv as args
from tqdm import tqdm
import re
from os.path import getsize

res = pd.DataFrame(columns=[
    "Window_Size",
    "N_Peaks",
    "Overlap",
    "Target_Zone",
    "Median_Hits",
    "Average_Hits",
    "Self_Matches",
    "Unique_Self_Matches",
    "Mean_F1_Score",
    "DB_Size_MB"
    ]
)
for file in tqdm(args[1:]):
    summary = pd.read_csv(file)
    grid_conf = list(re.findall('WINSIZE_(\d+)_NPEAKS_(\d+)_OVERLAP_(\d+)_DIFF_(\d+)',file)[0])
    grid_conf[-1] = 2 ** int(grid_conf[-1])
    res.loc[len(res.index)] = (
        *grid_conf,
        summary["First_Match_Count"].median(),
        round(summary["First_Match_Count"].mean(), 2),
        (summary["Sample_In_First_Matches"]).sum(),
        (summary["Sample_In_First_Matches"] & (summary["First_Match_Count"] == 1)).sum(),
        summary["F1_Score"].mean(),
        getsize("%s.pickle" % file.rsplit(".", 2)[0])/(1024**2)
    )

print(res.sort_values("Mean_F1_Score", ascending=False).to_csv(index=False, float_format="%g"))
