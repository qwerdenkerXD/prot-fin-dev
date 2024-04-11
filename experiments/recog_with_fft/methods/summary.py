import pandas as pd
from glob import glob
from sys import argv as args
import re

print("Filename,Median_Hits,Average_Hits,Self_Matches,Unique_Self_Matches")
for file in args[1:]:
    summary = pd.read_csv(file)
    print(
        file,
        summary["First_Match_Count"].median(),
        round(summary["First_Match_Count"].mean(), 2),
        (summary["Sample_In_First_Matches"]).sum(),
        (summary["Sample_In_First_Matches"] & (summary["First_Match_Count"] == 1)).sum(),
        sep=","
    )
