from .find_matches import get_filtered_db
from tools import *
import pandas as pd
from tqdm import tqdm


def match_family(
        family_file: str,
        db_in: str,
        filter_quantile=1.0
        ):

    (database, protein_lookup), hash_blacklist = get_filtered_db(db_in, filter_quantile)
    fam_data = pd.read_csv(family_file, sep=",", index_col="Protein_ID").squeeze()

    fam_hashes = {f: {} for f in fam_data.unique()}
    for hash_, prots in tqdm(database.items()):
        fam = fam_data.loc[[prot for _, prot in prots]]
        for f in fam:
            fam_hashes[f][hash_] = fam_hashes[f].get(hash_, 0) + 1

    print("Family_ID", "F_Score", "Precision", "Recall", "Member_Count", "Match_Count", "Hash_Intersec_Size", sep=",")
    for fam in tqdm(fam_hashes):
        fam_member_count = (fam_data == fam).sum()
        if fam_member_count > 1:
            hashes = set(hash_ for hash_, count in fam_hashes[fam].items() if count == fam_member_count)
            if len(hashes):
                match_prots = {}
                for hash_ in hashes:
                    for _, match_prot in database.get(hash_):
                        match_prots[match_prot] = match_prots.get(match_prot, 0) + 1

                match_fams = fam_data.loc[sorted(match_prots, key=lambda x: match_prots[x] / len(hashes))]
                true_pos_cumsum = (match_fams == fam).groupby(match_fams.index).max().cumsum()
                true_pos = true_pos_cumsum[-1]
                false_neg = fam_member_count - true_pos
                assert false_neg == 0
                positives = len(match_prots) - true_pos_cumsum.value_counts().get(0, 0)
                false_pos = positives - true_pos

                precision = true_pos / positives
                assert precision >= 0, (true_pos, positives)
                recall = true_pos / (true_pos + false_neg)

                f1_score = (2 * precision * recall) / (precision + recall)
                print(fam, round(f1_score, 2), round(precision, 2), round(recall, 2), fam_member_count, len(match_prots), len(hashes), sep=",")

            else:
                eprint("No intersection hashes for", fam)
        else:
            eprint("Ignored family", fam, "with only", fam_member_count, "member")
