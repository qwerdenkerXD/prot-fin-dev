from tools import *
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import pickle


def plot_family_covering(db: str, mapman: str, plot_out: str):
    with open(db, "rb") as f:
        database, _ = pickle.load(f)

    proteins = pd.read_csv(mapman, sep="\t", quotechar="'", index_col=1, usecols=["IDENTIFIER", "BINCODE"])
    proteins = proteins[proteins.index.notna()]
    proteins = proteins["BINCODE"].apply(str)

    family_covering = {i: {} for i in proteins.unique()}#
    for hash_, prots in tqdm(database.items()):
        fam = proteins.loc[[prot.lower() for _, prot in prots]]
        for f in fam:
            family_covering[f][hash_] = family_covering[f].get(hash_, 0) + 1

    plt.figure(figsize=(200, 5))

    print("Index", "Family", "Member_Count", "Max_Covering", "Max_Covering_Hashes", sep=",")
    for i, fam in enumerate(family_covering):
        fam_member_count = (proteins == fam).sum()
        counts = np.array(list(family_covering[fam].values()))
        plt.boxplot(counts / fam_member_count, positions=[i], widths=.8)
        max_cov = counts.max() if len(counts) else 0
        print(
            i,
            fam,
            fam_member_count,
            round(max_cov / fam_member_count, 2),
            ";".join(sorted(str(h) for h, c in family_covering[fam].items() if c == max_cov)),
            sep=","
        )

    plt.xticks(range(0, len(family_covering), 50), labels=range(0, len(family_covering), 50))
    ymin, ymax = plt.ylim()
    plt.ylim(ymin, 1 - ymin)
    plt.xlabel("Mapman Bin Index")
    plt.ylabel("Covering")
    plt.title("Covering of families by their members' hashes")
    plt.savefig(plot_out, bbox_inches='tight')
