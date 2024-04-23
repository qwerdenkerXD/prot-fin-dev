import pandas as pd
from sys import argv as args
import re

print("Filename,Median_Hits,Self_Matches,Unique_Self_Matches,Median_Score_Input_Family,Median_Score_Not_Input_Family,Mean_Score_Input_Family,Mean_Score_Not_Input_Family")
for file in args[1:]:
    summary = pd.read_csv(file)
    print(
        file,
        summary["First_Match_Count"].median(),
        (summary["Sample_In_First_Matches"]).sum(),
        (summary["Sample_In_First_Matches"] & (summary["First_Match_Count"] == 1)).sum(),
        round(summary["Median_Score_Input_Family"].mean(), 2),
        round(summary["Median_Score_Not_Input_Family"].mean(), 2),
        round(summary["Mean_Score_Input_Family"].mean(), 2),
        round(summary["Mean_Score_Not_Input_Family"].mean(), 2),
        sep=","
    )
