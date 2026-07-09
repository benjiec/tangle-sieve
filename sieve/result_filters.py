import importlib
import re


class RuleResultFilter(object):

    def __and__(self, other):
        return AndFilter(self, _as_filter(other))

    def __or__(self, other):
        return OrFilter(self, _as_filter(other))

    def __invert__(self):
        return NotFilter(self)

    def __call__(self, row):
        return self.matches(row)

    def matches(self, row):
        raise NotImplementedError

    def validate_columns(self, columns):
        raise NotImplementedError


def _as_filter(value):
    if not isinstance(value, RuleResultFilter):
        raise TypeError(f"Expected RuleResultFilter, got {type(value)}")
    return value


class CompositeFilter(RuleResultFilter):

    def __init__(self, *filters):
        self.filters = [_as_filter(filter_) for filter_ in filters]

    def validate_columns(self, columns):
        for filter_ in self.filters:
            filter_.validate_columns(columns)


class AndFilter(CompositeFilter):

    def matches(self, row):
        return all(filter_.matches(row) for filter_ in self.filters)


class OrFilter(CompositeFilter):

    def matches(self, row):
        return any(filter_.matches(row) for filter_ in self.filters)


class NotFilter(RuleResultFilter):

    def __init__(self, filter_):
        self.filter = _as_filter(filter_)

    def matches(self, row):
        return not self.filter.matches(row)

    def validate_columns(self, columns):
        self.filter.validate_columns(columns)


class FieldSelector(object):

    def eq(self, value):
        return FieldComparison(self, "eq", value)

    def ne(self, value):
        return FieldComparison(self, "ne", value)

    def num_eq(self, value):
        return FieldComparison(self, "num_eq", value)

    def num_ne(self, value):
        return FieldComparison(self, "num_ne", value)

    def gt(self, value):
        return FieldComparison(self, "gt", value)

    def gte(self, value):
        return FieldComparison(self, "gte", value)

    def ge(self, value):
        return self.gte(value)

    def lt(self, value):
        return FieldComparison(self, "lt", value)

    def lte(self, value):
        return FieldComparison(self, "lte", value)

    def le(self, value):
        return self.lte(value)

    def matches(self, pattern):
        return FieldComparison(self, "matches", pattern)

    def not_matches(self, pattern):
        return FieldComparison(self, "not_matches", pattern)

    def matching_columns(self, columns):
        raise NotImplementedError

    def compare_matches(self, row, operator, expected):
        columns = self.matching_columns(row.keys())
        if not columns:
            raise ValueError(f"Rule results are missing field matching: {self.describe()}")
        return any(_value_matches(row[column], operator, expected) for column in columns)

    def validate_columns(self, columns):
        if not self.matching_columns(columns):
            raise ValueError(f"Rule results are missing field matching: {self.describe()}")

    def describe(self):
        raise NotImplementedError


class Field(FieldSelector):

    def __init__(self, name):
        self.name = name

    def matching_columns(self, columns):
        return [column for column in columns if column == self.name]

    def describe(self):
        return self.name


class FieldRegex(FieldSelector):

    def __init__(self, pattern):
        self.pattern = pattern
        self.regex = re.compile(pattern)

    def matching_columns(self, columns):
        return [column for column in columns if self.regex.fullmatch(column) is not None]

    def describe(self):
        return self.pattern

    def all(self):
        return MultiFieldSelector(self, "all")

    def any(self):
        return MultiFieldSelector(self, "any")

    def compare_matches(self, row, operator, expected):
        raise ValueError(f"FieldRegex('{self.pattern}') requires .all() or .any() before comparison")


class MultiFieldSelector(FieldSelector):

    def __init__(self, selector, mode):
        self.selector = selector
        self.mode = mode

    def matching_columns(self, columns):
        return self.selector.matching_columns(columns)

    def compare_matches(self, row, operator, expected):
        columns = self.matching_columns(row.keys())
        if not columns:
            raise ValueError(f"Rule results are missing field matching: {self.describe()}")
        value_matches = [_value_matches(row[column], operator, expected) for column in columns]
        if self.mode == "all":
            return all(value_matches)
        if self.mode == "any":
            return any(value_matches)
        raise ValueError(f"Unknown multi-field mode: {self.mode}")

    def describe(self):
        return self.selector.describe()


def LeaderCall(prediction):
    return Field(f"Leader.call('{prediction}')")


class FieldComparison(RuleResultFilter):

    def __init__(self, selector, operator, expected):
        self.selector = selector
        self.operator = operator
        self.expected = expected

    def matches(self, row):
        return self.selector.compare_matches(row, self.operator, self.expected)

    def validate_columns(self, columns):
        self.selector.validate_columns(columns)


def _numeric(value, context):
    try:
        return float(value)
    except (TypeError, ValueError) as e:
        raise ValueError(f"Expected numeric value for {context}, got: {value!r}") from e


def _value_matches(actual, operator, expected):
    actual = "" if actual is None else str(actual)
    if operator == "eq":
        return actual == str(expected)
    if operator == "ne":
        return actual != str(expected)
    if operator == "matches":
        return re.search(str(expected), actual) is not None
    if operator == "not_matches":
        return re.search(str(expected), actual) is None
    actual_number = _numeric(actual, "rule row")
    expected_number = _numeric(expected, "result filter")
    if operator == "num_eq":
        return actual_number == expected_number
    if operator == "num_ne":
        return actual_number != expected_number
    if operator == "gt":
        return actual_number > expected_number
    if operator == "gte":
        return actual_number >= expected_number
    if operator == "lt":
        return actual_number < expected_number
    if operator == "lte":
        return actual_number <= expected_number
    raise ValueError(f"Unknown operator: {operator}")


def load_result_filter(spec):
    module_name, sep, attr_name = spec.rpartition(".")
    if not sep:
        raise ValueError("Result filter spec must be in the form module.attribute")

    module = importlib.import_module(module_name)
    value = getattr(module, attr_name)
    if isinstance(value, RuleResultFilter):
        return value
    raise TypeError(f"{spec} must refer to a sieve.result_filters.RuleResultFilter instance")
