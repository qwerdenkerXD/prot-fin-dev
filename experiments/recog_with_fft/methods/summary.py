import pandas as pd
from glob import glob
from sys import argv as args
from tqdm import tqdm
from tools import eprint
from os.path import exists, getsize
import re

res = pd.DataFrame(columns=[
    "Window_Size",
    "N_Peaks",
    "Overlap",
    "Hash_Use_First_Appearance",
    "Quantile",
    "Median_Hits_Top_Score",
    "Average_Hits",
    "Self_Matches",
    "Unique_Self_Matches",
    "Mean_F1_Score",
    "Mean_Precision",
    "Mean_Liberal_F1_Score",
    "Mean_Liberal_Precision",
    "Mean_Sharpness",
    "DB_Size_MB",
    "Runtime_find-matches",
    ]
)
for file in tqdm(args[1:]):
    try:
        summary = pd.read_csv(file)
        logfile = file.split("/")
        logfile.insert(-1, "_logs")
        logfile = glob("%s*.err" % "/".join(logfile).rsplit("_FILT", 1)[0])[0]
        params = re.findall('WINSIZE_(\d+)_NPEAKS_(\d+)_OVERLAP_(\d+)_FIRSTAPPEARANCE_(\d+)_FILT_(\.\d+)', file)[0]
        with open(logfile, "r") as f:
            log = f.readlines()
            time_regex = re.compile("([:\d]+)<[:\d]+,[ \d\.it/s\]]+$")
            for i, line in enumerate(log):
                if line == "findmatches FILT %s\n" % params[-1]:
                    db_size = re.search("\(([\d\.]+)MB\) of database", log[i+1]).group(1)
                elif line == "eval FILT %s\n" % params[-1]:
                    find_matches_time = time_regex.search(log[i-1]).group(1)
                    break  # eval comes after findmatches -> all conditions met
            del log

            if find_matches_time is not None:
                res.loc[len(res.index)] = (
                    *params,
                    summary["First_Match_Count"].median(),
                    round(summary["First_Match_Count"].mean(), 2),
                    (summary["Sample_In_First_Matches"]).sum(),
                    (summary["Sample_In_First_Matches"] & (summary["First_Match_Count"] == 1)).sum(),
                    summary["F1_Score"].mean(),
                    summary["Precision"].mean(),
                    summary["Liberal_F1_Score"].mean(),
                    summary["Precision_Liberal"].mean(),
                    summary["Sharpness"].mean(),
                    db_size,
                    find_matches_time,
                )
    except Exception as e:
        # raise e
        eprint("Problems with file '%s': May be empty: %s" % (file, e))

print(res.sort_values("Mean_F1_Score", ascending=False).to_csv(index=False, float_format="%g"))
