import importlib.util

from sieve.rules import RULE_FALSE, RULE_MAYBE, RULE_TRUE, Rule


def load_script(path):
    spec = importlib.util.spec_from_file_location("script_under_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ConstantByProteinRule(Rule):

    label = "ConstantByProteinRule"

    def evaluate(self, context):
        return {
            "p_true": RULE_TRUE,
            "p_maybe": RULE_MAYBE,
        }.get(context.protein.protein_accession, RULE_FALSE)


class AnnotatedByProteinRule(ConstantByProteinRule):

    label = "AnnotatedByProteinRule"

    def annotation_columns(self):
        return ["Example.call"]

    def annotations_many(self, contexts, rule_results):
        return {
            context.key: {"Example.call": context.protein.protein_accession + "_call"}
            for context in contexts
        }
