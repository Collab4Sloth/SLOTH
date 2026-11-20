#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Consistency checks of a SLOTH coefficient file (JSON), run before the generation of the C++
classes by GenerateCoefficient.py.

Errors (the generation is stopped):
    - structure: unknown or missing attributes, wrong types, two coefficients with the same
      class_name, both or none of 'expression' and 'expressions',
    - names: invalid identifiers, names declared twice (variables, explicit and auxiliary
      variables, coordinates, constants), C++ keywords and names used by the generated code, names
      read by SymPy as one of its objects (E, I, S, N, gamma, lambda...),
    - constants: format (name:value), value which is neither a number, nor R, K, NA, H, T, nor
      an arithmetic expression of numbers, of C++ constants M_* and of the other constants,
    - expressions: undeclared symbols (also after the expansion of the sums), unknown
      functions, dot()/sdot() outside a gradient coefficient, conditional expressions and
      bounds (same checks as the generation).

Warnings:
    - names of the arguments of the evaluation structure (values, time, F...): the coefficient
      is generated without this structure (no evaluation by batch nor on the device),
    - constants, explicit and auxiliary variables which are not used in the expressions (nor in
      the bounds). Variables and coordinates are not concerned: they may be declared without
      being used (e.g. a source term independent of the unknown).
"""
import difflib
import json
import keyword
import logging
import re
from collections import Counter
from dataclasses import dataclass
from typing import Dict, List

import sympy as sp
from sympy.core.function import AppliedUndef

from GenerateCoefficient import (EVAL_RESERVED_NAMES, expand_ranges, get_bounds,
                                 get_coordinates, read_coefficient)
from CoefficientReport import Dot, ExpressionParser, TIME_VALUE

# Attributes of a coefficient, and of a piece of a conditional expression
KNOWN_KEYS = {"class_name", "outputfile", "variables", "explicit_variables",
              "auxiliary_variables", "constants", "coordinates", "bounds", "gradient",
              "expression", "expressions", "range_variable"}
REQUIRED_KEYS = ("class_name", "outputfile", "variables")
STRING_KEYS = {"class_name", "outputfile", "variables", "explicit_variables",
               "auxiliary_variables", "constants", "coordinates", "range_variable", "expression"}
PIECE_KEYS = {"expression", "lower", "upper", "lower_strict", "upper_strict", "default"}

# Values of a constant other than a number (see get_constants)
CONSTANT_KEYWORDS = {"R", "K", "NA", "H", TIME_VALUE}

# Functions translated in C++
ALLOWED_FUNCTIONS = {"sin", "cos", "tan", "asin", "acos", "atan", "atan2", "sinh", "cosh", "tanh",
                     "asinh", "acosh", "atanh", "exp", "log", "sqrt", "Abs", "Max", "Min",
                     "sign", "floor", "ceiling", "erf", "erfc"}

# C++ keywords
CPP_KEYWORDS = {
    "alignas", "alignof", "and", "and_eq", "asm", "auto", "bitand", "bitor", "bool", "break",
    "case", "catch", "char", "class", "compl", "concept", "const", "consteval", "constexpr",
    "constinit", "const_cast", "continue", "co_await", "co_return", "co_yield", "decltype",
    "default", "delete", "do", "double", "dynamic_cast", "else", "enum", "explicit", "export",
    "extern", "false", "float", "for", "friend", "goto", "if", "inline", "int", "long",
    "mutable", "namespace", "new", "noexcept", "not", "not_eq", "nullptr", "operator", "or",
    "or_eq", "private", "protected", "public", "register", "reinterpret_cast", "requires",
    "return", "short", "signed", "sizeof", "static", "static_assert", "static_cast", "struct",
    "switch", "template", "this", "thread_local", "throw", "true", "try", "typedef", "typeid",
    "typename", "union", "unsigned", "using", "virtual", "void", "volatile", "wchar_t", "while",
    "xor", "xor_eq"}
# Names used by the generated C++ code: functions and namespaces called in the expressions, and
# arguments of the lambdas of the class written without evaluation structure
CPP_RESERVED = CPP_KEYWORDS | {"pow", "fmin", "fmax", "std", "mfem", "input_vector",
                               "exp_input_vector", "auxiliary_vector", "dimension"}

# Declared names reported if they are not used
WARN_IF_UNUSED = ("explicit variable", "auxiliary variable", "time", "constant")

IDENTIFIER = re.compile(r'^[A-Za-z_]\w*$')
CONSTANT_ITEM = re.compile(r'^\(\s*([^:()]*?)\s*:\s*([^:()]*?)\s*\)$')


@dataclass
class Issue:
    severity: str  # "error" or "warning"
    coefficient: str
    message: str

    def __str__(self):
        return f"{self.coefficient}: {self.message}"


class _Collector:
    def __init__(self, coefficient: str):
        self.coefficient = coefficient
        self.issues: List[Issue] = []

    def error(self, message: str):
        self.issues.append(Issue("error", self.coefficient, message))

    def warning(self, message: str):
        self.issues.append(Issue("warning", self.coefficient, message))

    @property
    def has_errors(self) -> bool:
        return any(i.severity == "error" for i in self.issues)


# ==========================================================
#  Checks
# ==========================================================

def _suggest(key: str, candidates) -> str:
    close = difflib.get_close_matches(key, sorted(candidates), n=1)
    return f" (did you mean '{close[0]}'?)" if close else ""


def _check_structure(item, out: _Collector):
    if not isinstance(item, dict):
        out.error("a coefficient must be a JSON object {...}")
        return
    for key in item:
        if key not in KNOWN_KEYS:
            out.error(f"unknown attribute '{key}'{_suggest(key, KNOWN_KEYS)}")
    for key in REQUIRED_KEYS:
        if key not in item:
            out.error(f"missing attribute '{key}'")
    if ("expression" in item) == ("expressions" in item):
        out.error("exactly one of 'expression' and 'expressions' must be defined")
    for key in STRING_KEYS & set(item):
        if not isinstance(item[key], str):
            out.error(f"'{key}' must be a string")
    if "gradient" in item and not isinstance(item["gradient"], bool):
        out.error("'gradient' must be true or false")
    if "bounds" in item and not isinstance(item["bounds"], dict):
        out.error("'bounds' must be an object {variable: [lower, upper]}")
    if "expressions" in item:
        pieces = item["expressions"]
        if not isinstance(pieces, list):
            out.error("'expressions' must be a list")
        else:
            for piece in pieces:
                if not isinstance(piece, dict):
                    out.error("each item of 'expressions' must be an object")
                    continue
                for key in piece:
                    if key not in PIECE_KEYS:
                        out.error(f"unknown attribute '{key}' in 'expressions'"
                                  f"{_suggest(key, PIECE_KEYS)}")


def _names(text: str) -> List[str]:
    expanded, _ = expand_ranges(text, 0)
    return [n for n in re.split(r'[,\s]+', expanded) if n]


def _parse_constants(text: str, out: _Collector) -> Dict[str, str]:
    constants = {}
    if not text.strip():
        return constants
    for token in re.findall(r'\([^()]*\)|[^,\s()][^,]*', text):
        match = CONSTANT_ITEM.match(token.strip())
        if not match or not match.group(1) or not match.group(2):
            out.error(f"constant '{token.strip()}' must be written (name:value)")
            continue
        name, value = match.groups()
        if name in constants:
            out.error(f"constant '{name}' defined twice")
        constants[name] = value
        if value not in CONSTANT_KEYWORDS and not _is_constant_value(value, constants):
            out.error(f"value '{value}' of the constant '{name}' must be a number, one of "
                      f"{', '.join(sorted(CONSTANT_KEYWORDS))}, or an arithmetic expression of "
                      "numbers, of C++ constants M_* and of the other constants")
    return constants


def _is_constant_value(value: str, constants: Dict[str, str]) -> bool:
    """
    Number, or arithmetic expression (copied in the C++ code) of numbers, of C++ constants M_*
    (M_PI...) and of the other constants, e.g. 8.*M_PI.
    """
    try:
        float(value)
        return True
    except ValueError:
        pass
    if not re.fullmatch(r'[\w.+\-*/() ]+', value):
        return False
    names = re.findall(r'[A-Za-z_]\w*', re.sub(r'\d+\.?\d*[eE][-+]?\d+', '', value))
    return all(n.startswith("M_") or n in constants for n in names)


def _sympy_meaning(name: str):
    """
    Object given by SymPy to a name read without declaration (as in the generation), or None
    if the name is read as a plain symbol.
    """
    try:
        obj = sp.sympify(name)
    except Exception:
        return "a reserved word of SymPy/Python"
    if isinstance(obj, sp.Symbol):
        return None
    known = {sp.E: "Euler's number e", sp.I: "the imaginary unit i", sp.pi: "the number pi"}
    if obj in known:
        return known[obj]
    return f"the SymPy object '{obj}' ({type(obj).__name__})"


def _check_names(item, out: _Collector) -> Dict[str, List[str]]:
    """
    Check the declared names and return them by category.
    """
    categories = {
        "variable": _names(item.get("variables", "")),
        "explicit variable": _names(item.get("explicit_variables", "")),
        "auxiliary variable": _names(item.get("auxiliary_variables", "")),
    }
    try:
        categories["coordinate"] = get_coordinates(item.get("coordinates", ""))
    except ValueError as e:
        out.error(str(e).replace("Error: ", "").replace(" Please check your data.", ""))
        categories["coordinate"] = []
    constants = _parse_constants(item.get("constants", ""), out)
    categories["time"] = [k for k, v in constants.items() if v == TIME_VALUE]
    categories["constant"] = [k for k, v in constants.items() if v != TIME_VALUE]

    if not categories["variable"]:
        out.error("at least one variable must be declared in 'variables'")

    counts = Counter(n for names in categories.values() for n in names)
    for name, count in counts.items():
        if count > 1:
            where = [c for c, names in categories.items() if name in names]
            out.error(f"'{name}' declared {count} times ({', '.join(where)})")

    for category, names in categories.items():
        for name in dict.fromkeys(names):
            if not IDENTIFIER.match(name):
                out.error(f"invalid {category} name '{name}'")
            elif name in CPP_RESERVED or keyword.iskeyword(name):
                out.error(f"{category} name '{name}' is reserved in the generated C++ code")
            elif name in EVAL_RESERVED_NAMES:
                out.warning(f"{category} name '{name}' is an argument of the evaluation "
                            "structure: the coefficient is generated without it (no evaluation "
                            "by batch nor on the device)")
            else:
                meaning = _sympy_meaning(name)
                if meaning is not None:
                    out.error(f"{category} name '{name}' is read by SymPy as {meaning}: "
                              "please rename it")

    class_name = item.get("class_name", "")
    if isinstance(class_name, str) and class_name and (
            not IDENTIFIER.match(class_name) or class_name in CPP_RESERVED):
        out.error(f"class_name '{class_name}' is not a valid C++ class name")
    return categories


def _check_expressions(item, spec, categories, out: _Collector) -> set:
    """
    Check the expressions and return the names they use.
    """
    declared = {n for names in categories.values() for n in names}
    used = set()
    parser = ExpressionParser(spec)
    for (expression, *_rest) in spec.expressions:
        if re.search(r"\bsum\s*\(", expression):
            out.error("sums are written Sum(...), not sum(...)")
            continue
        if not spec.is_gradient and re.search(r"\b(s?dot)\s*\(", expression):
            out.error("dot() and sdot() are only allowed in a gradient coefficient "
                      "(\"gradient\": true)")
            continue
        try:
            expr = parser.expanded(expression)
        except Exception as e:
            out.error(f"expression '{expression}' cannot be read: {e}")
            continue
        for f in expr.atoms(sp.Function):
            name = type(f).__name__
            if isinstance(f, Dot):
                continue
            if isinstance(f, AppliedUndef) or name not in ALLOWED_FUNCTIONS:
                out.error(f"unknown function '{name}' in '{expression}'")
        symbols = {str(s) for s in expr.free_symbols}
        for name in sorted(symbols - declared):
            out.error(f"'{name}' is used in '{expression}' but is not declared (check also the "
                      "bounds of the sums and the ranges a(i..j))")
        used |= symbols

    range_variable = item.get("range_variable")
    if range_variable is not None:
        if range_variable not in declared:
            out.error(f"range_variable '{range_variable}' is not declared")
        used.add(range_variable)
    return used


def _check_bounds(item, categories, out: _Collector) -> set:
    bounds = item.get("bounds")
    if not bounds:
        return set()
    if item.get("gradient", False):
        out.error("'bounds' are not allowed in a gradient coefficient")
        return set()
    variables = (categories["variable"] + categories["explicit variable"]
                 + categories["auxiliary variable"])
    try:
        get_bounds(bounds, variables, categories["constant"], out.coefficient)
    except ValueError as e:
        message = re.sub(r"^Error \([^)]*\): ", "", str(e)).replace(" Please check your data.", "")
        out.error(message)
    used = set()
    for value in bounds.values():
        if isinstance(value, (list, tuple)):
            for bound in value:
                if isinstance(bound, str):
                    used |= set(re.findall(r'[A-Za-z_]\w*', bound))
    return used


def check_coefficient(item, index: int) -> List[Issue]:
    name = item.get("class_name") if isinstance(item, dict) else None
    out = _Collector(name if isinstance(name, str) and name else f"coefficient #{index + 1}")

    _check_structure(item, out)
    if out.has_errors:
        return out.issues

    categories = _check_names(item, out)
    if out.has_errors:
        return out.issues

    try:
        spec = read_coefficient(item)
    except (ValueError, KeyError) as e:
        out.error(str(e).replace("Error: ", ""))
        return out.issues

    used = _check_expressions(item, spec, categories, out)
    used |= _check_bounds(item, categories, out)

    if not out.has_errors:
        for category in WARN_IF_UNUSED:
            for n in dict.fromkeys(categories[category]):
                if n not in used:
                    out.warning(f"{category} '{n}' is declared but not used")
    return out.issues


def check_file(input_file) -> List[Issue]:
    """
    Check a coefficient file and return the list of issues (errors and warnings).
    """
    try:
        with open(input_file, "r") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        return [Issue("error", str(input_file), f"invalid JSON: {e}")]
    except OSError as e:
        return [Issue("error", str(input_file), f"cannot be read: {e}")]

    if not isinstance(data, list):
        return [Issue("error", str(input_file), "the file must contain a list [...] of "
                                                "coefficients")]
    issues = []
    names = Counter(item.get("class_name") for item in data if isinstance(item, dict))
    for name, count in names.items():
        if name is not None and count > 1:
            issues.append(Issue("error", name, f"class_name defined {count} times in the file"))
    for index, item in enumerate(data):
        issues += check_coefficient(item, index)
    return issues


def log_issues(issues: List[Issue], input_file) -> bool:
    """
    Log the issues and return True if there is no error.
    """
    for issue in issues:
        if issue.severity == "error":
            logging.error(f"{input_file}: {issue}")
        else:
            logging.warning(f"{input_file}: {issue}")
    nb_errors = sum(i.severity == "error" for i in issues)
    if nb_errors:
        logging.error(f"{input_file}: {nb_errors} error(s), C++ files not generated.")
    return nb_errors == 0
