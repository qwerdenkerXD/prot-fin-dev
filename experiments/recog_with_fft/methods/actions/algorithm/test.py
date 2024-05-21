from test import TestCase
from .kidera import *
from .hash_gen import *
from .constellation import *
from tools import *


class TestKidera(TestCase):
    def test_extend_selected_factor(self):
        extended = self.create_valid(
            pd.Series,
            KIDERA_TABLE.iloc[0]
        )
        self.assertIsNone(extend_selected_factor(extended), "Unexpected return value")

        for special in "XBZJΨΩΦζΠ+-":
            self.assertIn(special, extended, "Missing special character in extended kidera table")

    def test_transform_seq(self):
        extended = self.create_valid(
            pd.Series,
            KIDERA_TABLE.iloc[0]
        )
        self.assertIsNone(extend_selected_factor(extended), "Unexpected return value")

        seq = "".join(extended.keys()) + "234567"  # numbers represent unknown amino acids
        transformed_seq = self.create_valid(
            List[np.float32],
            transform_seq(seq, extended, True)
        )
        for aa, val in zip(seq, transformed_seq):
            self.assertEqual(val, np.float32(extended.get(aa, np.float32(0))), "Sequence falsey transformed")


class TestConstellation(TestCase):
    def test_find_peaks(self):
        spectrum = [10, 9, 0, 1, 0, 7, 6]
        peaks = [5, 3]
        returned = self.create_valid(
            List[int],
            find_peaks(spectrum, 1)
        )
        self.assertEqual(returned, peaks[:1])
        returned = self.create_valid(
            List[int],
            find_peaks(spectrum, 0)
        )
        self.assertEqual(returned, peaks)


class TestHashGen(TestCase):
    def test_create_hash(self):
        hash_ = self.create_valid(
            Hash,
            create_hash((10, 10), (20, 20))
        )
        self.assertEqual(hash_ % 2**20, 20, "Falsey hash creation")
        self.assertEqual(hash_ >> 20, 10, "Falsey hash creation")

        self.assertRaisesRegex(AssertionError, "exceeds 32 bit", create_hash, (0, 33))
        self.assertRaisesRegex(AssertionError, "too big for 1 bit", create_hash, (3, 1))

        const_map = self.create_valid(
            ConstellationMap,
            [((0, 0., 0),), ((0, 0., 0),)]
        )
        hashes = self.create_valid(
            Hashes,
            create_hashes(const_map, ProteinID(), 0)
        )
        self.assertEqual(len(hashes), 2, "Expected exactly two created hash")
        expected_hashes = (
            Hash(0),  # diff is 0, freqs are 0, factor is 0, quantiles are 0 -> 0
            Hash(2 ** FREQUENCY_BITS - 1 << FREQUENCY_BITS),  # last frequency will be combined with dummy, as there are no upcoming ones
        )
        for i, expected_hash in enumerate(expected_hashes):
            self.assertIn(expected_hash, hashes, "Wrong hash")
            self.assertEqual(hashes[expected_hash], (i, ProteinID()))
