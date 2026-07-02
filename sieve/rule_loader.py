import importlib

from sieve.rules import Rule, Rules


def load_rules(spec):
    module_name, sep, attr_name = spec.rpartition(".")
    if not sep:
        raise ValueError("Rule spec must be in the form module.attribute")

    module = importlib.import_module(module_name)
    value = getattr(module, attr_name)
    if isinstance(value, Rules):
        return value
    if isinstance(value, Rule):
        return Rules(value)
    raise TypeError(f"{spec} must refer to a sieve.rules.Rules or Rule instance")
