import os
import sys
import tempfile
import unittest

from sieve.rule_loader import load_rules
from sieve.rules import Rules


class TestRuleLoader(unittest.TestCase):

    def test_load_rules_accepts_rules_or_rule_instances(self):
        with tempfile.TemporaryDirectory() as tmpd:
            module_path = os.path.join(tmpd, "example_rules.py")
            with open(module_path, "w", encoding="utf-8") as f:
                f.write("\n".join([
                    "from sieve.rules import KO, Rules",
                    "raw_rule = KO.matches('K04564')",
                    "wrapped_rule = Rules(raw_rule)",
                    "",
                ]))
            sys.path.insert(0, tmpd)
            try:
                self.assertIsInstance(load_rules("example_rules.raw_rule"), Rules)
                self.assertIsInstance(load_rules("example_rules.wrapped_rule"), Rules)
            finally:
                sys.path.remove(tmpd)
                sys.modules.pop("example_rules", None)
