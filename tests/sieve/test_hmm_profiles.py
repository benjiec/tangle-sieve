import os
import tempfile
import unittest

from sieve.hmm_profiles import hmm_profiles_in_dirs
from sieve.rules import _hmm_profiles_by_basename


class TestHMMProfiles(unittest.TestCase):

    def test_finds_only_nonrecursive_hmm_files_in_sorted_order(self):
        with tempfile.TemporaryDirectory() as tmpd:
            nested = os.path.join(tmpd, "nested")
            os.makedirs(nested)
            for path in [
                os.path.join(tmpd, "b.hmm"),
                os.path.join(tmpd, "a.hmm"),
                os.path.join(tmpd, "ignored.txt"),
                os.path.join(nested, "nested.hmm"),
            ]:
                with open(path, "w", encoding="utf-8"):
                    pass

            self.assertEqual(
                hmm_profiles_in_dirs([tmpd]),
                [os.path.join(tmpd, "a.hmm"), os.path.join(tmpd, "b.hmm")],
            )

    def test_rejects_missing_directory(self):
        with tempfile.TemporaryDirectory() as tmpd:
            missing = os.path.join(tmpd, "missing")
            with self.assertRaisesRegex(ValueError, "HMM profile directory does not exist"):
                hmm_profiles_in_dirs([missing])

    def test_duplicate_basenames_are_rejected_by_profile_registry(self):
        with tempfile.TemporaryDirectory() as tmpd:
            first = os.path.join(tmpd, "first")
            second = os.path.join(tmpd, "second")
            os.makedirs(first)
            os.makedirs(second)
            for directory in [first, second]:
                with open(os.path.join(directory, "profile.hmm"), "w", encoding="utf-8"):
                    pass

            profiles = hmm_profiles_in_dirs([first, second])
            with self.assertRaisesRegex(ValueError, "Ambiguous HMM profile basename"):
                _hmm_profiles_by_basename(profiles)
