from test import TestCase
from tools import *
import pandas as pd
import pickle
import sys
import os

from .create_db import *
from .find_matches import *


class TestCreateDB(TestCase):
    protein_file = "test/create_db.fa"
    db_out = "test/create_db.pickle"

    @classmethod
    def setUpClass(cls):
        cls.assertIsNone(cls, create_db(cls.protein_file, cls.db_out))

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.db_out)

    def test_create_db(self):
        with open(self.db_out, "rb") as f:
            db = self.create_valid(
                DB,
                pickle.load(f)
            )

        self.assertEqual(len(db.lookup), 3, "Protein lookup is not complete")

        for hash_, idx_prot_pairs in db.db.items():

            for pair in idx_prot_pairs:
                self.assertEqual(len(pair), 2, "Protein index pair in DB has different length")
                _, prot = pair
                self.assertIn(prot, db.lookup, f"'{prot}' not in protein lookup")

        for prot_id, prot_info in db.lookup.items():
            self.assertEqual(len(prot_info), 2, "Lookup value is differs in length")

            description, hash_count = prot_info
            self.assertGreaterEqual(hash_count, 0, "Hash count in protein lookup below zero")


class TestFindMatches(TestCase):
    protein_file = "test/create_db.fa"
    db_in = "test/find_matches.pickle"
    stdout_pipe = "test/find_matches.matches.tmp"

    @classmethod
    def setUpClass(cls):
        cls.assertIsNone(cls, create_db(cls.protein_file, cls.db_in))
        with open(cls.db_in, "rb") as f:
            cls.db, cls.lookup = pickle.load(f)[:2]

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.db_in)

    def tearDown(self):
        sys.stdout = sys.__stdout__
        if os.path.exists(self.stdout_pipe):
            os.remove(self.stdout_pipe)

    def test_find_matches(self):
        # Just testing if the function runs without errors

        with open(self.stdout_pipe, "w") as f:
            sys.stdout = f
            self.assertIsNone(find_matches(self.protein_file, self.db_in))

    def test_score_prots(self):
        with open(self.db_in, "rb") as f:
            db, lookup = pickle.load(f)[:2]

        hashes: Hashes = self.create_valid(
            Hashes,
            {hash_: (WindowIndex(i), ProteinID()) for i, hash_ in enumerate(db.keys())}
        )
        scores: ScoresMap = self.create_valid(
            ScoresMap,
            score_prots(hashes, db, lookup)
        )
        self.assertEqual(db, self.db, "Database has changed")
        self.assertEqual(lookup, self.lookup, "Database has changed")

    def test_get_matches_per_prot(self):
        hashes: Hashes = self.create_valid(
            Hashes,
            {Hash(123): (WindowIndex(1), ProteinID())}
        )
        proteins = [(WindowIndex(0), ProteinID(i)) for i in range(10)]
        db = self.create_valid(
            Database,
            {Hash(123): proteins}
        )
        matches: MatchesPerProt = self.create_valid(
            MatchesPerProt,
            get_matches_per_prot(hashes, db)
        )
        self.assertEqual(len(matches), len(proteins), "Returned matches not complete")
        for idx, prot_id in proteins:
            self.assertIn(prot_id, matches, f"Protein identifier '{prot_id}' missing")
            original_offsets = self.create_valid(
                ScoresByOffset,
                {WindowIndex(1): 1}
            )
            self.assertEqual(
                matches[prot_id],
                original_offsets,
                f"Found offsets for protein id '{prot_id}' wrong"
            )

    def test_get_max_offset(self):
        scores = list(range(5)) + list(range(3, 0, -1))
        scored_offsets: ScoresByOffset = self.create_valid(
            ScoresByOffset,
            {WindowIndex(i): score for i, score in enumerate(scores)}
        )
        max_offset: (WindowIndex, Score) = self.create_valid(
            Tuple[WindowIndex, Score],
            get_max_offset(scored_offsets)
        )
        max_score = sorted(scores)[-1]
        self.assertEqual(max_offset, (WindowIndex(scores.index(max_score)), max_score))

    def test_get_result_frame(self):
        proteins = [ProteinID(i) for i in range(10)]
        scores_map: ScoresMap = self.create_valid(
            ScoresMap,
            {p: (WindowIndex(i), Score(10 - + i//2), JSI(1 - i//3 * .1)) for i, p in enumerate(proteins)}
        )
        result: pd.DataFrame = self.create_valid(
            pd.DataFrame,
            get_result_frame(scores_map, proteins[0], 10, 5)
        )
        self.assertTrue(get_result_frame({}, proteins[0], 10, 5).empty, "Expected empty dataframe for zero scores")

        rank = 0
        prev_score, prev_jsi = float("+inf"), float("+inf")

        for expected, actual in zip(scores_map.items(), result.sort_values("Rank").values):
            if actual[2] < prev_jsi or actual[3] < prev_score:
                prev_jsi, prev_score = actual[2], actual[3]
                rank += 1
            (e_match_id, (_, e_score, e_jsi)) = expected
            e_input_id, e_seq_len, e_hash_count = proteins[0], 10, 5
            self.assertEqual(
                (rank, e_match_id, e_jsi, e_score, proteins[0], 10, 5),
                tuple(actual),
                "Protfin output not correct"
            )


class TestEvaluateProtfin(TestCase):
    ...


class TestSelectSamples(TestCase):
    ...
