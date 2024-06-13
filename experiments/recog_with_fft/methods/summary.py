import pandas as pd
from glob import glob
from sys import argv as args
from tools import eprint
import re

res = pd.DataFrame(columns=[
    "Window_Size",
    "N_Peaks",
    "Overlap",
    "Quantile",
    "Median_Hits",
    "Average_Hits",
    "Self_Matches",
    "Unique_Self_Matches",
    "Mean_F1_Score",
    "DB_Size_MB"
    ]
)
for file in args[1:]:
    try:
        summary = pd.read_csv(file)
        logfile = file.split("/")
        logfile.insert(-1, "_logs")
        logfile = glob("%s*.err" % "/".join(logfile).rsplit(".", 2)[0])[0]
        with open(logfile, "r") as f:
            res.loc[len(res.index)] = (
                *re.findall('WINSIZE_(\d+)_NPEAKS_(\d+)_OVERLAP_(\d+)_FILT_(\.\d+)', file)[0],
                summary["First_Match_Count"].median(),
                round(summary["First_Match_Count"].mean(), 2),
                (summary["Sample_In_First_Matches"]).sum(),
                (summary["Sample_In_First_Matches"] & (summary["First_Match_Count"] == 1)).sum(),
                summary["F1_Score"].mean(),
                re.findall("([\d.]+MB)\) of database size", f.read())[0]
            )
    except Exception as e:
        eprint("Problems with file '%s': May be empty" % file)

print(res.sort_values("Mean_F1_Score", ascending=False).to_csv(index=False, float_format="%g"))
