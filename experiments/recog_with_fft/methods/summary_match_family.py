import pandas as pd
from glob import glob
from sys import argv as args
from tqdm import tqdm
from tools import eprint
import numpy as np
import re

res = pd.DataFrame(columns=[
    "Window_Size",
    "N_Peaks",
    "Overlap",
    "Hash_Use_First_Appearance",
    "Quantile",
    "Average_F_Score",
    "Average_Precision",
    "avg(Family_Members/Match_Count)",
    "Runtime_match-family"
    ]
)
for file in tqdm(args[1:]):
    try:
        summary = pd.read_csv(file)
        logfile = file.split("/")
        logfile.insert(-1, "_logs")
        logfile = glob("%s*.err" % "/".join(logfile).rsplit("_FILT", 1)[0])[0]
        params = re.findall('WINSIZE_(\d+)_NPEAKS_(\d+)_OVERLAP_(\d+)_FIRSTAPPEARANCE_(\d+)_FILT_(\.?\d+)', file)[0]
        f_scores = summary["F_Score"][~summary["F_Score"].isna()]

        with open(logfile, "r") as f:
            log = f.readlines()
            create_db_time = None
            for i, line in enumerate(log):
                if line == "matchfamily FILT %s\n" % params[-1]:
                    match_family_time = re.search("([:\d]+)<[:\d]+,[ \d\.it/s\]]+$", log[i-1])
                    break
            del log

            if match_family_time is not None:
                res.loc[len(res.index)] = (
                    *params,
                    np.pad(f_scores, (0, 4937-len(f_scores))).mean(),
                    summary["Precision"].mean(),
                    (summary["Member_Count"] / summary["Match_Count"]).mean(),
                    match_family_time.group(1)
                )
    except Exception as e:
        # raise e
        eprint("Problems with file '%s': May be empty" % file)

res["Hash_Use_First_Appearance"] = res["Hash_Use_First_Appearance"].apply(lambda x: bool(int(x)))
print(res.sort_values("Average_F_Score", ascending=False).to_csv(index=False, float_format="%g"))
