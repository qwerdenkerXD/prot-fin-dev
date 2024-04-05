import pandas as pd
from glob import glob
import re

print("Window_Size,Overlap,Selected_Peaks,Median_Hits,Average_Hits,Self_Matches,Unique_Self_Matches")
for file in sorted(glob("../results/_stft_param_exp/*.summary.csv")):
    summary = pd.read_csv(file)
    print(
          *re.findall(".*_(\d+).*_(\d+).*_(\d+)", file)[0],
          summary["First_Match_Count"].median(),
          round(summary["First_Match_Count"].mean(), 2),
          (summary["Sample_In_First_Matches"]).sum(),
          (summary["Sample_In_First_Matches"] & (summary["First_Match_Count"] == 1)).sum(),
          sep=","
          )
