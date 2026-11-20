#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Variational coefficients of SLOTH: energy densities depending on the gradients of the variables.

A variational coefficient ("variational": true in the JSON file) defines an energy density
f(u, grad(u)), u being the variables of the coefficient. The components of the gradients are
handled as independent symbols p_{i,k} = du_i/dx_k, so that the first variation of F = int f
in the direction v, i.e. the weak form of the Euler-Lagrange equation, reads
    dF[u; v] = int ( df/du_i v_i + df/dp_{i,k} dv_i/dx_k )
and its second variation (Jacobian) involves d2f/du2, d2f/dudp and d2f/dp2.

Syntax of the expressions:
    grad(a)         gradient of a (variable, explicit or auxiliary variable)
    n(a)            normal grad(a) / sqrt(|grad(a)|^2 + eps^2), eps given by "regularization"
    dot(A, B)       scalar product of A and B, each one being grad(.) or n(.)
    A[k]            component k (0 <= k < dimension) of A = grad(.) or n(.)
    sdot(a(i..j))   dot(grad(ai), grad(ai)) + ... + dot(grad(aj), grad(aj))

Gradients of the variables are differentiated. Gradients of the explicit and auxiliary
variables are data: they are not differentiated.

Storage: component k of the gradient of the variable i has the index i*dimension + k.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import sympy as sp

# Default value of eps in n(a) = grad(a) / sqrt(|grad(a)|^2 + eps^2)
DEFAULT_REGULARIZATION = 1.0e-12

# Prefix of the symbols standing for the components of the gradients: grad_<name>_<k>
GRADIENT_PREFIX = "grad_"

# Names of the components k = 0, 1, 2 in the LaTeX report
AXES = ("x", "y", "z")

# Names of the functions available in the expressions of a variational coefficient
VARIATIONAL_FUNCTIONS = ("grad", "n", "dot", "sdot")


class VariationalError(ValueError):
    """
    Error in the expression of a variational coefficient.
    """


# ==========================================================
#  Symbolic objects read in the expressions
# ==========================================================

def axis(k) -> str:
    """
    Name of the component k (x, y, z) in the LaTeX report.
    """
    integer = isinstance(k, int) or getattr(k, "is_Integer", False)
    return AXES[int(k)] if integer and 0 <= int(k) < len(AXES) else str(k)


def _power(tex: str, exp) -> str:
    """
    LaTeX of tex raised to the power exp (given by the printer for a power of a function).
    """
    return tex if exp is None else f"\\left({tex}\\right)^{{{exp}}}"


class Grad(sp.Function):
    """
    grad(a): gradient of a.
    """
    nargs = 1

    def __getitem__(self, k):
        return Component(self, k)

    def _latex(self, printer, exp=None):
        return _power(f"\\nabla {printer._print(self.args[0])}", exp)


class Normal(sp.Function):
    """
    n(a): regularized unit normal grad(a) / sqrt(|grad(a)|^2 + eps^2).
    """
    nargs = 1

    def __getitem__(self, k):
        return Component(self, k)

    def _latex(self, printer, exp=None):
        return _power(f"\\mathbf{{n}}_{{{printer._print(self.args[0])}}}", exp)


class Component(sp.Function):
    """
    A[k]: component k of A = grad(a) or n(a).
    """
    nargs = 2

    def _latex(self, printer, exp=None):
        vector, k = self.args
        name, k = printer._print(vector.args[0]), axis(k)
        if isinstance(vector, Normal):
            tex = f"\\mathbf{{n}}_{{{name},{k}}}"
            return tex if exp is None else f"{tex}^{{{exp}}}"
        return _power(f"\\partial_{{{k}}} {name}", exp)


class VDot(sp.Function):
    """
    dot(A, B): scalar product of A and B, each one being grad(a) or n(a).
    """
    nargs = 2

    def _latex(self, printer, exp=None):
        a, b = [printer._print(arg) for arg in self.args]
        if self.args[0] == self.args[1]:
            return _power(f"\\left|{a}\\right|^{{2}}", exp)
        return _power(f"{a} \\cdot {b}", exp)


# Local names given to the parser of the expressions of a variational coefficient
PARSER_FUNCTIONS = {"grad": Grad, "n": Normal, "dot": VDot}

VARIATIONAL_CLASSES = (Grad, Normal, Component, VDot)


def gradient_symbols(name: str, dimension: int) -> List[sp.Symbol]:
    """
    Symbols of the components of the gradient of name: grad_<name>_0, ..., grad_<name>_<d-1>.
    """
    return [sp.Symbol(f"{GRADIENT_PREFIX}{name}_{k}") for k in range(dimension)]


# ==========================================================
#  Expression in terms of the components of the gradients
# ==========================================================

def to_components(expr, dimension: int, regularization: float, gradient_names: Sequence[str]):
    """
    Replace grad(.), n(.), dot(., .) and the components A[k] in expr by expressions of the
    symbols grad_<name>_<k>.

    regularization: eps (number, or positive SymPy symbol),
    gradient_names: names whose gradient can be taken (variables, explicit and auxiliary
    variables).
    Return (expression, names whose gradient is used).
    """
    eps = regularization if isinstance(regularization, sp.Basic) else sp.Float(regularization)
    eps2 = eps ** 2
    allowed = set(gradient_names)
    used = set()

    def vector(obj):
        if not isinstance(obj, (Grad, Normal)):
            raise VariationalError(f"'{obj}' is neither grad(.) nor n(.): dot() and [k] apply "
                                   "to grad(a) or n(a) only")
        arg = obj.args[0]
        if not isinstance(arg, sp.Symbol):
            raise VariationalError(f"'{obj}': the argument of grad() and n() must be a name")
        if arg.name not in allowed:
            raise VariationalError(f"'{obj}': '{arg.name}' is not a variable, an explicit or an "
                                   "auxiliary variable")
        used.add(arg.name)
        g = gradient_symbols(arg.name, dimension)
        if isinstance(obj, Grad):
            return g
        norm = sp.sqrt(sum(c ** 2 for c in g) + eps2)
        return [c / norm for c in g]

    def dot(a, b):
        return sp.Add(*[x * y for x, y in zip(vector(a), vector(b))])

    def component(obj, k):
        if not (k.is_Integer and 0 <= int(k) < dimension):
            raise VariationalError(f"'{obj}[{k}]': the component must be an integer between 0 "
                                   f"and {dimension - 1} (dimension {dimension})")
        return vector(obj)[int(k)]

    expr = expr.replace(lambda e: isinstance(e, VDot), lambda e: dot(*e.args))
    expr = expr.replace(lambda e: isinstance(e, Component), lambda e: component(*e.args))
    leftover = expr.atoms(Grad, Normal)
    if leftover:
        names = ", ".join(sorted(str(e) for e in leftover))
        raise VariationalError(f"{names}: grad(.) and n(.) must be used in dot(.,.) or with a "
                               "component [k]")
    return expr, used


# ==========================================================
#  Derivatives
# ==========================================================

@dataclass
class VariationalForm:
    """
    Energy density f and its derivatives with respect to the variables u_i and to the
    components p_{i,k} of their gradients (index i*dimension + k).

        df_dvar          [n]           df/du_i
        df_dgrad         [n*d]         df/dp_{i,k}
        d2f_dvar2        [n x n]       d2f/du_i du_j
        d2f_dvar_dgrad   [n x n*d]     d2f/du_i dp_{j,k}
        d2f_dgrad2       [n*d x n*d]   d2f/dp_{i,k} dp_{j,l}
    """
    variables: List[sp.Symbol]
    dimension: int
    gradients: List[sp.Symbol]
    data_gradients: Dict[str, List[sp.Symbol]]
    density: sp.Expr
    df_dvar: List[sp.Expr] = field(default_factory=list)
    df_dgrad: List[sp.Expr] = field(default_factory=list)
    d2f_dvar2: sp.Matrix = None
    d2f_dvar_dgrad: sp.Matrix = None
    d2f_dgrad2: sp.Matrix = None

    @property
    def nb_variables(self) -> int:
        return len(self.variables)

    def depends_on_gradients(self) -> bool:
        """
        True if f depends on the gradient of at least one variable.
        """
        return any(self.density.has(p) for p in self.gradients)

    def isotropic_coefficients(self) -> Optional[List[sp.Expr]]:
        """
        Coefficients lambda_i if the gradient part of f is isotropic and decoupled:
            df/dp_{i,k} = lambda_i p_{i,k},  d2f/dp2 = diag(lambda_i I_d),  d2f/du dp = 0,
        lambda_i depending on no gradient (f = sum_i lambda_i/2 |grad(u_i)|^2 + g(u)).
        Return None otherwise.
        """
        d = self.dimension
        all_gradients = set(self.gradients)
        for g in self.data_gradients.values():
            all_gradients |= set(g)
        if any(sp.simplify(e) != 0 for e in self.d2f_dvar_dgrad):
            return None
        lambdas = []
        for i in range(self.nb_variables):
            lam = sp.simplify(self.d2f_dgrad2[i * d, i * d])
            if lam.free_symbols & all_gradients:
                return None
            lambdas.append(lam)
        size = self.nb_variables * d
        for r in range(size):
            for c in range(size):
                expected = lambdas[r // d] if r == c else 0
                if sp.simplify(self.d2f_dgrad2[r, c] - expected) != 0:
                    return None
            if sp.simplify(self.df_dgrad[r] - lambdas[r // d] * self.gradients[r]) != 0:
                return None
        return lambdas


def build_variational_form(expr, variable_names: Sequence[str], dimension: int,
                           regularization: float = DEFAULT_REGULARIZATION,
                           data_names: Sequence[str] = ()) -> VariationalForm:
    """
    Variational form of the energy density expr (SymPy expression with grad, n, dot, A[k]).

    variable_names: variables (differentiated, with their gradients),
    regularization: eps of n(.) (number, or positive SymPy symbol),
    data_names: explicit and auxiliary variables (their gradients are data).
    """
    if dimension not in (1, 2, 3):
        raise VariationalError(f"dimension must be 1, 2 or 3, not {dimension}")
    positive = regularization.is_positive if isinstance(regularization, sp.Basic) \
        else regularization > 0
    if not positive:
        raise VariationalError(f"regularization must be positive, not {regularization}")
    density, used = to_components(expr, dimension, regularization,
                                  list(variable_names) + list(data_names))
    variables = [sp.Symbol(v) for v in variable_names]
    gradients = [p for v in variable_names for p in gradient_symbols(v, dimension)]
    data_gradients = {v: gradient_symbols(v, dimension) for v in data_names if v in used}

    form = VariationalForm(variables, dimension, gradients, data_gradients, density)
    form.df_dvar = [sp.diff(density, u) for u in variables]
    form.df_dgrad = [sp.diff(density, p) for p in gradients]
    form.d2f_dvar2 = sp.Matrix([[sp.diff(fu, v) for v in variables] for fu in form.df_dvar])
    form.d2f_dvar_dgrad = sp.Matrix([[sp.diff(fu, p) for p in gradients]
                                     for fu in form.df_dvar])
    form.d2f_dgrad2 = sp.Matrix([[sp.diff(fp, q) for q in gradients] for fp in form.df_dgrad])
    return form


def variational_form_of(spec, expression: Optional[str] = None) -> VariationalForm:
    """
    Variational form of a coefficient read from the JSON file (CoefficientSpec).
    expression: expression to use (default: the expression of the coefficient).
    """
    from CoefficientReport import ExpressionParser
    parser = ExpressionParser(spec)
    if expression is None:
        expression = spec.expressions[0][0]
    return build_variational_form(parser.expanded(expression), parser.variables,
                                  spec.dimension, spec.regularization,
                                  parser.explicit_variables + parser.auxiliary_variables)
