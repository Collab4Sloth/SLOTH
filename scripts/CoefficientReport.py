#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LaTeX report describing the analytical expressions of a SLOTH coefficient file (JSON).

The report is built from the CoefficientSpec objects read by GenerateCoefficient.py, so that
it describes exactly the data used to generate the C++ classes. It is a description of the
JSON file: it does not modify the generated C++ code.

For each coefficient (class_name), the report gives:
    - the kind of coefficient (standard, conditional, gradient),
    - the variables (differentiated), the explicit and auxiliary variables, the coordinates,
    - the constants and their values,
    - the bounds of the variables,
    - the expression F, with sums Sum(...) kept in compact form,
    - the first derivatives of F and its second derivatives (if the Hessian is small enough).
"""
import logging
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import sympy as sp
from sympy.parsing.sympy_parser import convert_xor, standard_transformations
from sympy.printing.latex import LatexPrinter

from VariationalForm import (PARSER_FUNCTIONS, Component, Normal, VDot, axis,
                             build_variational_form, gradient_symbols)

RANGE_PATTERN = re.compile(r'^([a-zA-Z_]\w*)\((\d+)\.\.(\d+)\)$')
SDOT_PATTERN = re.compile(r"sdot\(\s*(\w+)\((\d+)\.\.(\d+)\)\s*\)")
INDEXED_PATTERN = re.compile(r'([a-zA-Z_]\w*)\s*\[')

# Hessian written in full only for at most this number of variables
MAX_HESSIAN_SIZE = 3
# Maximal visible width (in characters, LaTeX commands not counted) of a displayed line:
# a longer sum is split over several lines, each line being filled with as many terms as
# possible
MAX_LINE_WIDTH = 60

TRANSFORMATIONS = standard_transformations + (convert_xor,)

PHYSICAL_CONSTANTS = {
    "R": "ideal gas constant (\\texttt{Physical::R})",
    "K": "Boltzmann constant (\\texttt{Physical::K})",
    "NA": "Avogadro constant (\\texttt{Physical::NA})",
    "H": "Planck constant (\\texttt{Physical::H})",
}

# Value of a constant standing for the current time (this->time_ in the C++ code)
TIME_VALUE = "T"


# ==========================================================
#  Symbols and parsing
# ==========================================================

class Dot(sp.Function):
    """
    dot(a, b) of the gradient coefficients: product of the gradients of a and b.
    """
    nargs = 2

    def _latex(self, printer):
        a, b = [printer._print(arg) for arg in self.args]
        if self.args[0] == self.args[1]:
            return f"\\left|\\nabla {a}\\right|^{{2}}"
        return f"\\nabla {a} \\cdot \\nabla {b}"


class ReportPrinter(LatexPrinter):
    """
    LaTeX printer with an explicit thin space between factors (names such as xe, ye or dvol
    are not mistaken for products of single letters).
    """

    def __init__(self, settings=None):
        settings = dict(settings or {})
        settings.setdefault("order", "none")
        settings.setdefault("mul_symbol", r"\,")
        super().__init__(settings)

    def _needs_mul_brackets(self, expr, first=False, last=False) -> bool:
        if isinstance(expr, (Dot, VDot)):
            return False
        return super()._needs_mul_brackets(expr, first=first, last=last)


def to_latex(expr) -> str:
    return ReportPrinter().doprint(expr)


def split_names(names: str) -> List[str]:
    """
    Split a list of names separated by commas and/or spaces, without expanding the ranges.
    """
    return [t for t in re.split(r'[,\s]+', names.strip()) if t]


def expand_names(names: str) -> List[str]:
    """
    Expand a(i..j) in ai, ..., aj.
    """
    result = []
    for token in split_names(names):
        match = RANGE_PATTERN.match(token)
        if match:
            name, start, end = match.groups()
            result += [f"{name}{i}" for i in range(int(start), int(end) + 1)]
        else:
            result.append(token)
    return result


def names_to_latex(names: str) -> str:
    """
    LaTeX of a list of names: a(1..9) is written a_1, ..., a_9.
    """
    items = []
    for token in split_names(names):
        match = RANGE_PATTERN.match(token)
        if match:
            name, start, end = match.groups()
            base = to_latex(sp.Symbol(name))
            start, end = int(start), int(end)
            if end - start <= 1:
                items += [f"{base}_{{{i}}}" for i in range(start, end + 1)]
            else:
                items.append(f"{base}_{{{start}}}, \\dots, {base}_{{{end}}}")
        else:
            items.append(to_latex(sp.Symbol(token)))
    return ", ".join(f"${i}$" for i in items)


def parse_constants(constants: str) -> Dict[str, str]:
    """
    Parse "(a:1.0),(R:R)" as an ordered dict name -> value (string).
    """
    if not constants:
        return {}
    result = {}
    for t in constants.split(','):
        t = t.replace("(", "").replace(")", "").strip()
        if not t:
            continue
        name, value = t.split(':')
        result[name.strip()] = value.strip()
    return result


class ExpressionParser:
    """
    Parse the expressions of one coefficient with every declared name as a plain symbol (no
    clash with SymPy names such as E, S, N, gamma...).
    """

    def __init__(self, spec):
        self.spec = spec
        self.variables = expand_names(spec.variables)
        self.explicit_variables = expand_names(spec.explicit_variables)
        self.auxiliary_variables = expand_names(spec.auxiliary_variables)
        self.coordinates = split_names(spec.coordinates.replace(",", " ")) \
            if spec.coordinates else []
        constants = parse_constants(spec.constants)
        # Time: constants of value T (current time), the other ones are true constants
        self.time = [k for k, v in constants.items() if v == TIME_VALUE]
        self.constants = {k: v for k, v in constants.items() if v != TIME_VALUE}

        names = (self.variables + self.explicit_variables + self.auxiliary_variables
                 + self.coordinates + self.time + list(self.constants))
        self.symbols = {n: sp.Symbol(n) for n in names}

    def _locals(self, expression: str) -> dict:
        local_dict = dict(self.symbols)
        if self.spec.is_gradient:
            local_dict["dot"] = Dot
        if self.spec.is_variational:
            local_dict.update(PARSER_FUNCTIONS)
        if "Sum(" in expression:
            local_dict["Sum"] = sp.Sum
            for base in set(INDEXED_PATTERN.findall(expression)):
                local_dict[base] = sp.IndexedBase(base)
            for index in set(re.findall(r'\[\s*([a-zA-Z_]\w*)\s*\]', expression)):
                local_dict[index] = sp.Symbol(index, integer=True)
        return local_dict

    def _prepare(self, expression: str) -> str:
        # sdot(eta(1..9)) -> Sum(dot(eta[i_],eta[i_]),(i_,1,9)) for a gradient coefficient,
        #                    Sum(dot(grad(eta[i_]),grad(eta[i_])),(i_,1,9)) for a variational one
        variational = self.spec.is_variational

        def repl(match):
            var, start, end = match.groups()
            term = f"grad({var}[i_])" if variational else f"{var}[i_]"
            return f"Sum(dot({term},{term}),(i_,{start},{end}))"
        return SDOT_PATTERN.sub(repl, expression)

    def display(self, expression: str):
        """
        Expression of the JSON file, with sums kept in compact form.
        """
        expression = self._prepare(expression)
        return sp.parse_expr(expression, local_dict=self._locals(expression),
                             transformations=TRANSFORMATIONS)

    def expanded(self, expression: str):
        """
        Expression with sums expanded on the named variables (eta[1] -> eta1), as in the
        generated C++ code. Used for the derivatives.
        """
        expr = self.display(expression).doit()
        indexed = {i: sp.Symbol(f"{i.base}{i.indices[0]}") for i in expr.atoms(sp.Indexed)}
        return expr.xreplace(indexed)

    def bound_to_latex(self, bound) -> Optional[str]:
        if bound is None:
            return None
        if isinstance(bound, (int, float)):
            return to_latex(sp.Float(bound)) if isinstance(bound, float) else str(bound)
        return to_latex(sp.parse_expr(bound, local_dict=dict(self.symbols), evaluate=False,
                                      transformations=TRANSFORMATIONS))


# ==========================================================
#  LaTeX writing
# ==========================================================

def escape(text: str) -> str:
    return re.sub(r'([_&%#$])', r'\\\1', str(text))


FUNCTION_NAMES = re.compile(r'\\(cos|sin|tan|exp|log|ln|sinh|cosh|tanh|sqrt)\b')
SPACING_COMMANDS = re.compile(r'\\(left|right|,|;|quad|qquad|operatorname)')


def visible_width(tex: str) -> int:
    """
    Rough width of a LaTeX formula: function names count for their length, other commands
    (greek letters, \\sum...) for one character, braces and spacing are not counted.
    """
    tex = SPACING_COMMANDS.sub("", tex)
    tex = FUNCTION_NAMES.sub(lambda m: m.group(1), tex)
    tex = re.sub(r'\\[a-zA-Z]+', "x", tex)
    tex = re.sub(r'[{}^_\s]', "", tex)
    return len(tex)


def split_expression(expr) -> List[str]:
    r"""
    LaTeX of expr split in lines of visible width up to MAX_LINE_WIDTH (as far as possible):
        - sum: lines filled with terms; a term too long alone is itself split,
        - product: split inside its longest factor which is a sum, with the delimiters
          \left( ... \right. and \left. ... \right) around the lines.
    Other expressions are not split.
    """
    tex = to_latex(expr)
    if visible_width(tex) <= MAX_LINE_WIDTH:
        return [tex]
    if isinstance(expr, sp.Add):
        return split_sum(expr)
    if isinstance(expr, sp.Mul):
        return split_product(expr)
    return [tex]


def split_sum(expr) -> List[str]:
    lines, current, width = [], "", 0
    for k, t in enumerate(expr.args):
        term = to_latex(t)
        negative = term.lstrip().startswith("-")
        sign = "" if (k == 0 or negative) else "+ "
        if visible_width(term) > MAX_LINE_WIDTH:
            # term too long alone: on its own lines
            if current:
                lines.append(current)
                current, width = "", 0
            sub = split_expression(t)
            lines += [sign + sub[0]] + sub[1:]
            continue
        piece = f"{sign}{term}" if k == 0 else f" {sign}{term}"
        w = visible_width(piece)
        if current and width + w > MAX_LINE_WIDTH:
            lines.append(current)
            current, width = "", 0
            piece = piece.lstrip()
        current += piece
        width += w
    if current:
        lines.append(current)
    return [line.strip() for line in lines]


def split_product(expr) -> List[str]:
    coeff, rest = expr.as_coeff_Mul()
    factors = sp.Mul.make_args(rest)
    sums = [f for f in factors if isinstance(f, sp.Add)]
    if not sums:
        return split_factors(coeff, factors)
    longest = max(sums, key=lambda f: visible_width(to_latex(f)))
    inner = split_expression(longest)
    if len(inner) == 1:
        return split_factors(coeff, factors)
    others = [f for f in factors if f is not longest]
    prefix = to_latex(coeff * sp.Mul(*others, evaluate=False)) if (others or coeff != 1) else ""
    if prefix == "1":
        prefix = ""
    elif prefix == "-1":
        prefix = "-"
    if prefix and prefix != "-":
        prefix += r"\,"
    lines = [f"{prefix}\\left({inner[0]}\\right."]
    lines += [f"\\left.{{}} {line}\\right." for line in inner[1:-1]]
    lines.append(f"\\left.{{}} {inner[-1]}\\right)")
    return lines


def split_factors(coeff, factors) -> List[str]:
    """
    Product whose factors cannot be split: lines filled with whole factors, a continuation
    line starting with \\times.
    """
    texs = [f"\\left({to_latex(f)}\\right)" if isinstance(f, sp.Add) else to_latex(f)
            for f in factors]
    sign = ""
    if coeff == -1:
        sign = "-"
    elif coeff != 1:
        texs.insert(0, to_latex(coeff))
    lines, current = [], sign
    for tex in texs:
        piece = tex if current in ("", "-") else f"\\,{tex}"
        if current not in ("", "-") and visible_width(current + piece) > MAX_LINE_WIDTH:
            lines.append(current)
            current, piece = "", f"\\times {tex}"
        current += piece
    lines.append(current)
    return lines


def equation(lhs: str, expr) -> List[str]:
    """
    Displayed equation lhs = expr (eqnarray). A long expression is split over several lines
    (see split_expression).
    """
    lines = split_expression(expr)
    rows = [f"  {lhs} &=& {lines[0]}"]
    rows += [f"  && \\quad {{}} {line}" for line in lines[1:]]
    return ["\\begin{eqnarray*}", " \\\\\n".join(rows), "\\end{eqnarray*}"]


def condition_to_latex(var: str, lower, lower_strict, upper, upper_strict, parser) -> str:
    v = to_latex(sp.Symbol(var))
    left = f"{parser.bound_to_latex(lower)} \\le " if lower is not None \
        else f"{parser.bound_to_latex(lower_strict)} < "
    right = f" \\le {parser.bound_to_latex(upper)}" if upper is not None \
        else f" < {parser.bound_to_latex(upper_strict)}"
    return f"{left}{v}{right}"


def write_expression(spec, parser) -> List[str]:
    out = []
    lhs = "F"
    if not spec.is_conditional:
        out += equation(lhs, parser.display(spec.expressions[0][0]))
        return out

    pieces = []
    default = None
    for (expr, var, lower, lower_strict, upper, upper_strict) in spec.expressions:
        tex = to_latex(parser.display(expr))
        if lower is None and lower_strict is None and upper is None and upper_strict is None:
            default = tex
        else:
            cond = condition_to_latex(var, lower, lower_strict, upper, upper_strict, parser)
            pieces.append(f"    {tex} & \\text{{if }} {cond}")
    pieces.append(f"    {default} & \\text{{otherwise}}")
    out += ["\\begin{eqnarray*}", f"  {lhs} &=&", "  \\begin{cases}",
            " \\\\\n".join(pieces), "  \\end{cases}", "\\end{eqnarray*}"]
    return out


def write_derivatives(spec, parser) -> List[str]:
    if spec.is_gradient:
        return ["Derivatives are not computed for gradient coefficients."]
    if spec.is_variational:
        return write_variational_derivatives(spec, parser)
    variables = [sp.Symbol(v) for v in parser.variables]
    n = len(variables)
    pieces = [parser.expanded(expr) for (expr, *_rest) in spec.expressions]

    def piece_title(k):
        return [f"\\paragraph{{Piece {k + 1}.}}"] if spec.is_conditional else []

    first = ["\\subsubsection{First derivatives}", ""]
    for k, f in enumerate(pieces):
        first += piece_title(k)
        for v in variables:
            first += equation(f"\\frac{{\\partial F}}{{\\partial {to_latex(v)}}}",
                              sp.diff(f, v))

    second = ["", "\\subsubsection{Second derivatives}", ""]
    if n > MAX_HESSIAN_SIZE:
        second.append(f"Hessian: ${n} \\times {n}$ matrix (symmetric), not written.")
        return first + second
    for k, f in enumerate(pieces):
        second += piece_title(k)
        # Hessian: upper triangle (symmetric), one entry per equation
        for i, vi in enumerate(variables):
            for vj in variables[i:]:
                if vi == vj:
                    lhs = f"\\frac{{\\partial^2 F}}{{\\partial {to_latex(vi)}^2}}"
                else:
                    lhs = (f"\\frac{{\\partial^2 F}}{{\\partial {to_latex(vi)}"
                           f"\\,\\partial {to_latex(vj)}}}")
                second += equation(lhs, sp.diff(f, vi, vj))
    return first + second


# ==========================================================
#  Variational coefficients
# ==========================================================

# Variational form written in full (second derivatives) for at most this number of unknowns
# (variables and components of their gradients)
MAX_VARIATIONAL_SIZE = 4
# Expressions are factorized for display up to this number of operations
MAX_TIDY_OPERATIONS = 2000


class Partial(sp.Function):
    """
    Display of the component k of the gradient of a: partial_x a.
    """
    nargs = 2

    def _latex(self, printer, exp=None):
        name, k = self.args
        tex = f"\\partial_{{{axis(k)}}} {printer._print(name)}"
        return tex if exp is None else f"\\left({tex}\\right)^{{{exp}}}"


class Modulus(sp.Function):
    """
    Display of the norm of the gradient of a: |grad a|.
    """
    nargs = 1

    def _latex(self, printer, exp=None):
        tex = f"\\left|\\nabla {printer._print(self.args[0])}\\right|"
        return tex if exp is None else f"{tex}^{{{exp}}}"


class VariationalDisplay:
    """
    Writing of the derivatives of a variational coefficient in the notations of the JSON file:
        - gradient of a name used without n(): components partial_x a, partial_y a...
        - gradient of a name used with n(): |grad a| n_{a,x}, in the limit eps -> 0 of the
          regularization (n = grad a / |grad a|).
    """

    def __init__(self, names, normal_names, dimension: int, eps):
        self.eps = eps
        self.sums = []
        self.replacements = {}
        for name in names:
            symbol = sp.Symbol(name)
            components = gradient_symbols(name, dimension)
            if name in normal_names:
                modulus = Modulus(symbol)
                self.sums.append((sum(c ** 2 for c in components), modulus ** 2))
                for k, c in enumerate(components):
                    self.replacements[c] = modulus * Component(Normal(symbol), k)
            else:
                for k, c in enumerate(components):
                    self.replacements[c] = Partial(symbol, k)

    def __call__(self, expr):
        expr = sp.sympify(expr).subs(self.eps, 0)
        for total, modulus2 in self.sums:
            expr = expr.subs(total, modulus2)
        expr = expr.xreplace(self.replacements)
        if sp.count_ops(expr) <= MAX_TIDY_OPERATIONS:
            # decimal numbers of the JSON file (0.5...) as fractions before factorization
            expr = sp.factor(sp.expand(sp.nsimplify(expr, rational=True)))
        return expr


def _partial(name: str, k=None) -> str:
    """
    LaTeX of the derivation variable: u, or the component k of grad u.
    """
    u = to_latex(sp.Symbol(name))
    return u if k is None else f"(\\partial_{{{axis(k)}}} {u})"


def write_variational_derivatives(spec, parser) -> List[str]:
    """
    Weak form of the variational coefficient: derivatives of f with respect to the variables
    and to the components of their gradients, strong form, natural boundary condition and
    second derivatives (Jacobian).
    """
    eps = sp.Dummy("epsilon", positive=True)
    d = spec.dimension
    names = parser.variables
    data = parser.explicit_variables + parser.auxiliary_variables
    expr = parser.expanded(spec.expressions[0][0])
    normal_names = {str(e.args[0]) for e in expr.atoms(Normal)}
    form = build_variational_form(expr, names, d, eps, data)
    display = VariationalDisplay(names + list(form.data_gradients), normal_names, d, eps)

    out = ["The energy density $f$ defines $F = \\int_\\Omega f \\, d\\Omega$. The weak form of "
           "the Euler-Lagrange equations is the first variation of $F$ in the direction $\\psi$:"]
    out += ["\\begin{equation*}",
            "  dF[\\boldsymbol{u};\\boldsymbol{\\psi}] = \\int_\\Omega \\sum_i \\left("
            "\\frac{\\partial f}{\\partial u_i}\\, \\psi_i + \\frac{\\partial f}{\\partial "
            "\\nabla u_i} \\cdot \\nabla \\psi_i \\right) d\\Omega",
            "\\end{equation*}"]
    out.append(f"with $u_i$ = {names_to_latex(spec.variables)}. Strong form (Euler-Lagrange) "
               "and natural boundary condition, $\\boldsymbol{\\nu}$ being the outward unit normal "
               "to the boundary:")
    out += ["\\begin{equation*}",
            "  \\frac{\\delta F}{\\delta u_i} = \\frac{\\partial f}{\\partial u_i} - \\nabla "
            "\\cdot \\frac{\\partial f}{\\partial \\nabla u_i}, \\qquad \\frac{\\partial f}"
            "{\\partial \\nabla u_i} \\cdot \\boldsymbol{\\nu} = 0 \\ \\text{on } "
            "\\partial\\Omega.",
            "\\end{equation*}",
            "A boundary term of the weak form prescribes the conormal flux $\\partial f / \\partial "
            "\\nabla u_i \\cdot \\boldsymbol{\\nu}$ (e.g. $g$ for a Neumann condition of value $g$), "
            "and not $\\nabla u_i \\cdot \\boldsymbol{\\nu}$: for an isotropic gradient term "
            "$\\partial f / \\partial \\nabla u_i = \\lambda_i \\nabla u_i$, it is "
            "$\\lambda_i \\nabla u_i \\cdot \\boldsymbol{\\nu}$."]
    if form.data_gradients:
        out.append(" The gradients of " + ", ".join(
            f"${to_latex(sp.Symbol(v))}$" for v in form.data_gradients)
            + " are data: they are not differentiated.")
    if normal_names:
        out.append(" Derivatives are written in the limit $\\varepsilon \\to 0$ of the "
                   "regularization of the normal, $\\mathbf{n}_a = \\nabla a / |\\nabla a|$; "
                   "the generated code uses $\\mathbf{n}_a = \\nabla a / \\sqrt{|\\nabla a|^2 + "
                   f"\\varepsilon^2}}$ with $\\varepsilon = {to_latex(sp.Float(spec.regularization))}$.")
    out.append("")

    lambdas = form.isotropic_coefficients()
    for i, name in enumerate(names):
        u = to_latex(sp.Symbol(name))
        out += [f"\\subsubsection{{Variable ${u}$}}", ""]
        out += equation(f"\\frac{{\\partial f}}{{\\partial {u}}}", display(form.df_dvar[i]))
        for k in range(d):
            out += equation(f"\\frac{{\\partial f}}{{\\partial {_partial(name, k)}}}",
                            display(form.df_dgrad[i * d + k]))
        if lambdas is not None:
            out.append(f"Isotropic gradient term: $\\partial f / \\partial \\nabla {u} = "
                       f"\\lambda \\nabla {u}$ with $\\lambda = {to_latex(display(lambdas[i]))}$.")
        out.append("")

    out += write_variational_hessian(form, names, display)
    return out


def write_variational_hessian(form, names, display) -> List[str]:
    """
    Second derivatives of f (Jacobian): non-zero entries of the upper triangle of the Hessian
    with respect to the variables and to the components of their gradients.
    """
    d = form.dimension
    unknowns = [(name, None) for name in names] + [(name, k) for name in names
                                                     for k in range(d)]
    size = len(unknowns)
    out = ["\\subsubsection{Second derivatives (Jacobian)}", ""]
    if size > MAX_VARIATIONAL_SIZE:
        out.append(f"Hessian with respect to the variables and to the components of their "
                   f"gradients: ${size} \\times {size}$ matrix (symmetric), not written.")
        return out
    n = len(names)

    def entry(r, c):
        if r < n and c < n:
            return form.d2f_dvar2[r, c]
        if r < n:
            return form.d2f_dvar_dgrad[r, c - n]
        return form.d2f_dgrad2[r - n, c - n]

    zeros = []
    for r in range(size):
        for c in range(r, size):
            value = display(entry(r, c))
            a, b = _partial(*unknowns[r]), _partial(*unknowns[c])
            lhs = (f"\\frac{{\\partial^2 f}}{{\\partial {a}^2}}" if r == c else
                   f"\\frac{{\\partial^2 f}}{{\\partial {a}\\,\\partial {b}}}")
            if value == 0:
                zeros.append(f"${lhs}$")
            else:
                out += equation(lhs, value)
    if zeros:
        out.append("Zero second derivatives: " + ", ".join(zeros) + ".")
    return out


def write_coefficient(spec) -> List[str]:
    parser = ExpressionParser(spec)
    out = [f"\\section{{Coefficient \\texttt{{{escape(spec.class_name)}}}}}", ""]

    kind = "conditional (piecewise)" if spec.is_conditional else \
        "gradient (post-processing only)" if spec.is_gradient else \
        f"variational (dimension {spec.dimension})" if spec.is_variational \
        else "standard"
    rows = [("Kind", kind),
            ("C++ class", f"\\texttt{{{escape(spec.class_name)}}} in "
                          f"\\texttt{{{escape(spec.output_file)}.hpp}}")]
    out += ["\\begin{tabular}{@{}ll@{}}"]
    out += [f"  \\textbf{{{k}}} & {v} \\\\" for k, v in rows]
    out += ["\\end{tabular}", ""]

    variables = [("Primary variables", names_to_latex(spec.variables))]
    if spec.explicit_variables:
        variables.append(("Explicit variables", names_to_latex(spec.explicit_variables)))
    if spec.auxiliary_variables:
        variables.append(("Auxiliary variables", names_to_latex(spec.auxiliary_variables)))
    if parser.coordinates:
        variables.append(("Coordinates", ", ".join(
            f"${to_latex(sp.Symbol(c))} = x_{{{i + 1}}}$"
            for i, c in enumerate(parser.coordinates))))
    if parser.time:
        variables.append(("Time", ", ".join(
            f"${to_latex(sp.Symbol(t))}$" if t == "t" else f"${to_latex(sp.Symbol(t))} = t$"
            for t in parser.time)))
    out += ["\\subsection{Variables}", "",
            "\\begin{center}",
            "\\begin{tabular}{@{}ll@{}}", "  \\toprule", "  Type & Symbols \\\\",
            "  \\midrule"]
    out += [f"  {k} & {v} \\\\" for k, v in variables]
    out += ["  \\bottomrule", "\\end{tabular}", "\\end{center}", ""]

    if parser.constants:
        out += ["\\subsection{Constants}", "",
                "\\begin{center}",
                "\\begin{tabular}{@{}ll@{}}", "  \\toprule", "  Symbol & Value \\\\",
                "  \\midrule"]
        for name, value in parser.constants.items():
            text = PHYSICAL_CONSTANTS.get(value, f"${to_latex(sp.Float(value))}$"
                                          if re.fullmatch(r'[-+.\deE]+', value) else
                                          f"\\texttt{{{escape(value)}}}")
            out.append(f"  ${to_latex(sp.Symbol(name))}$ & {text} \\\\")
        out += ["  \\bottomrule", "\\end{tabular}", 
                "\\end{center}",""]

    if spec.bounds:
        out += ["\\subsection{Bounds}", "",
                "Variables are clamped to their bounds before the evaluation of $F$ and of its "
                "derivatives:", "", "\\begin{itemize}"]
        for key, (lower, upper) in spec.bounds.items():
            v = names_to_latex(key)
            lo, up = parser.bound_to_latex(lower), parser.bound_to_latex(upper)
            if lo is not None and up is not None:
                text = f"${lo} \\le$ {v} $\\le {up}$"
            elif lo is not None:
                text = f"{v} $\\ge {lo}$"
            else:
                text = f"{v} $\\le {up}$"
            out.append(f"  \\item {text}")
        out += ["\\end{itemize}", ""]

    out += ["\\subsection{Expression}", ""]
    out += write_expression(spec, parser)
    out.append("")

    out += ["\\subsection{Derivatives with respect to the variables}", ""]
    out += write_derivatives(spec, parser)
    out.append("")
    return out


def write_report(coefficients, tex_file, json_file) -> Path:
    """
    Write the LaTeX report of the coefficients read from json_file.
    """
    tex_file = Path(tex_file).resolve()
    out = [
        "\\documentclass[11pt]{article}",
        "\\usepackage[margin=2cm]{geometry}",
        "\\usepackage{amsmath,amssymb,booktabs}",
        "\\allowdisplaybreaks",
        "\\setlength{\\parindent}{0pt}",
        "\\title{SLOTH coefficients: \\texttt{" + escape(Path(json_file).name) + "}}",
        "\\date{}",
        "\\begin{document}",
        "\\maketitle",
        "",
        "This report describes the analytical expressions defined in \\texttt{"
        + escape(Path(json_file).name) + "}. Each coefficient $F$ depends on the variables "
        "(unknowns, with respect to which $F$ is differentiated), on the explicit and "
        "auxiliary variables, on the spatial coordinates, on the time and on constants.",
        "",
    ]
    for spec in coefficients:
        out += write_coefficient(spec)
    out.append("\\end{document}")
    tex_file.write_text("\n".join(out) + "\n")
    logging.info(f"Report written in {tex_file}")
    return tex_file


def compile_report(tex_file: Path) -> Optional[Path]:
    """
    Compile the report with pdflatex if available. Return the PDF path or None.
    """
    pdflatex = shutil.which("pdflatex")
    if pdflatex is None:
        logging.warning("pdflatex not found: report not compiled.")
        return None
    cmd = [pdflatex, "-interaction=nonstopmode", "-halt-on-error", tex_file.name]
    result = subprocess.run(cmd, cwd=tex_file.parent, capture_output=True, text=True)
    stem = tex_file.with_suffix("")
    if result.returncode != 0:
        logging.warning(f"pdflatex failed, see {stem}.log")
        return None
    for ext in (".aux", ".log"):
        Path(f"{stem}{ext}").unlink(missing_ok=True)
    pdf = tex_file.with_suffix(".pdf")
    logging.info(f"Report compiled in {pdf}")
    return pdf
