#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests of the variational coefficients (scripts/VariationalForm.py) and of their checks
(scripts/CoefficientChecker.py).

Run from the build directory with ctest (ctest -L Scripts), or directly:
    PYTHONPATH=<SLOTH>/scripts python3 -m unittest -v test_VariationalForm
"""
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import sympy as sp

from CoefficientChecker import check_coefficient, check_file
from GenerateCoefficient import read_coefficient
from VariationalForm import (DEFAULT_REGULARIZATION, VariationalError, build_variational_form,
                             gradient_symbols, variational_form_of)

SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "scripts"
TESTS_DIR = Path(__file__).resolve().parents[1]

# Cubic (four-fold) anisotropy of PFHub benchmark 3, written with the components of n(phi)
ANISOTROPY = "(1 + e4*(n(phi)[0]**4 - 6*n(phi)[0]**2*n(phi)[1]**2 + n(phi)[1]**4))"


def coefficient(expression, variables="phi", dimension=2, **attributes):
    item = {"class_name": "Test", "outputfile": "Test", "variables": variables,
            "variational": True, "dimension": dimension, "expression": expression}
    item.update(attributes)
    return item


def form(expression, **kwargs):
    return variational_form_of(read_coefficient(coefficient(expression, **kwargs)))


def errors(item):
    return [i.message for i in check_coefficient(item, 0) if i.severity == "error"]


def warnings(item):
    return [i.message for i in check_coefficient(item, 0) if i.severity == "warning"]


def is_zero(expr) -> bool:
    if isinstance(expr, sp.MatrixBase):
        return all(sp.simplify(e) == 0 for e in expr)
    return sp.simplify(expr) == 0


class TestDerivatives(unittest.TestCase):
    """
    Derivatives of the energy density with respect to the variables and to the components of
    their gradients.
    """

    def test_isotropic(self):
        f = form("0.5*kappa*dot(grad(phi), grad(phi))", constants="(kappa:2.0)")
        kappa = sp.Symbol("kappa")
        p = gradient_symbols("phi", 2)
        self.assertEqual(f.gradients, p)
        self.assertTrue(all(is_zero(e) for e in f.df_dvar))
        self.assertTrue(all(is_zero(f.df_dgrad[k] - kappa * p[k]) for k in range(2)))
        self.assertTrue(is_zero(f.d2f_dgrad2 - kappa * sp.eye(2)))
        self.assertTrue(f.d2f_dvar_dgrad.is_zero_matrix)
        lambdas = f.isotropic_coefficients()
        self.assertEqual(len(lambdas), 1)
        self.assertTrue(is_zero(lambdas[0] - kappa))

    def test_sizes_follow_dimension(self):
        for d in (1, 2, 3):
            f = form("dot(grad(a), grad(b)) + a*b", variables="a, b", dimension=d)
            self.assertEqual(len(f.df_dvar), 2)
            self.assertEqual(len(f.df_dgrad), 2 * d)
            self.assertEqual(f.d2f_dvar2.shape, (2, 2))
            self.assertEqual(f.d2f_dvar_dgrad.shape, (2, 2 * d))
            self.assertEqual(f.d2f_dgrad2.shape, (2 * d, 2 * d))

    def test_coefficient_depending_on_the_variable(self):
        # f = (k0 + k1 phi)/2 |grad(phi)|^2: df/dphi = k1/2 |grad(phi)|^2, d2f/dphi dp = k1 p
        f = form("0.5*(k0 + k1*phi)*dot(grad(phi), grad(phi))", constants="(k0:1.0),(k1:0.5)")
        k1 = sp.Symbol("k1")
        p = gradient_symbols("phi", 2)
        self.assertTrue(is_zero(f.df_dvar[0] - k1 * (p[0] ** 2 + p[1] ** 2) / 2))
        self.assertTrue(all(is_zero(f.d2f_dvar_dgrad[0, k] - k1 * p[k]) for k in range(2)))
        self.assertIsNone(f.isotropic_coefficients())

    def test_cross_gradient_terms(self):
        # f = 1/2 (k11 |grad a|^2 + 2 k12 grad a . grad b + k22 |grad b|^2): d2f/dp2 = K x I_d
        f = form("0.5*(k11*dot(grad(a), grad(a)) + 2*k12*dot(grad(a), grad(b))"
                 " + k22*dot(grad(b), grad(b)))", variables="a, b", dimension=3,
                 constants="(k11:2.0),(k12:0.5),(k22:3.0)")
        K = sp.Matrix(sp.symbols("k11 k12 k12 k22")).reshape(2, 2)
        expected = sp.kronecker_product(K, sp.eye(3))
        self.assertTrue(is_zero(f.d2f_dgrad2 - expected))
        self.assertIsNone(f.isotropic_coefficients())

    def test_sdot(self):
        f = form("k*sdot(eta(1..3))/2", variables="eta(1..3)", constants="(k:2.0)")
        k = sp.Symbol("k")
        self.assertEqual(len(f.gradients), 6)
        self.assertTrue(is_zero(f.d2f_dgrad2 - k * sp.eye(6)))
        self.assertEqual(f.isotropic_coefficients(), [k, k, k])

    def test_sum_of_gradients(self):
        f = form("Sum(dot(grad(eta[i]), grad(eta[i])), (i, 1, 2))", variables="eta(1..2)")
        self.assertTrue(is_zero(f.d2f_dgrad2 - 2 * sp.eye(4)))

    def test_gradient_of_an_auxiliary_variable_is_data(self):
        # Coupling term lambda(c)/2 |grad(phi)|^2 in the equation of c (phi auxiliary):
        # df/dc = l/2 |grad(phi)|^2, no derivative with respect to grad(phi)
        f = form("0.5*l*c*dot(grad(phi), grad(phi))", variables="c",
                 auxiliary_variables="phi", constants="(l:2.0)")
        q = gradient_symbols("phi", 2)
        self.assertEqual(f.gradients, gradient_symbols("c", 2))
        self.assertEqual(f.data_gradients, {"phi": q})
        self.assertTrue(is_zero(f.df_dvar[0] - sp.Symbol("l") * (q[0] ** 2 + q[1] ** 2) / 2))
        self.assertTrue(all(is_zero(e) for e in f.df_dgrad))
        self.assertFalse(f.depends_on_gradients())
        self.assertEqual(f.isotropic_coefficients(), [0])

    def test_normal_is_regularized(self):
        f = form("n(phi)[0]", regularization=1.0e-3)
        p = gradient_symbols("phi", 2)
        expected = p[0] / sp.sqrt(p[0] ** 2 + p[1] ** 2 + sp.Float(1.0e-3) ** 2)
        self.assertTrue(is_zero(f.density - expected))
        self.assertEqual(read_coefficient(coefficient("dot(grad(phi), grad(phi))")).regularization,
                         DEFAULT_REGULARIZATION)


class TestAnisotropy(unittest.TestCase):
    """
    Four-fold anisotropy of PFHub benchmark 3: f = 1/2 W(n)^2 |grad(phi)|^2,
    W(n) = W0 (1 + e4 cos(4 theta)) with cos(4 theta) = nx^4 - 6 nx^2 ny^2 + ny^4.
    """

    def anisotropic_form(self, e4):
        return form(f"0.5*(W0*{ANISOTROPY})**2*dot(grad(phi), grad(phi))",
                    constants=f"(W0:1.0),(e4:{e4})")

    @staticmethod
    def values(e4, theta, modulus=1.0):
        p = gradient_symbols("phi", 2)
        return {sp.Symbol("W0"): 1.0, sp.Symbol("e4"): e4,
                p[0]: modulus * math.cos(theta), p[1]: modulus * math.sin(theta)}

    def test_flux_matches_the_benchmark(self):
        # df/dp = W^2 p + |p|^2 W dW/dp, with the exact normal p/|p|
        f = self.anisotropic_form(0.05)
        px, py = sp.symbols("px py")
        modulus = sp.sqrt(px ** 2 + py ** 2)
        nx, ny = px / modulus, py / modulus
        W = 1 + sp.Rational(1, 20) * (nx ** 4 - 6 * nx ** 2 * ny ** 2 + ny ** 4)
        reference = [sp.diff(W ** 2 * modulus ** 2 / 2, v) for v in (px, py)]
        values = self.values(0.05, 0.7, modulus=0.8)
        p = gradient_symbols("phi", 2)
        point = {px: values[p[0]], py: values[p[1]]}
        for k in range(2):
            computed = float(f.df_dgrad[k].subs(values))
            expected = float(reference[k].subs(point))
            self.assertAlmostEqual(computed, expected, places=12)

    def test_symmetry(self):
        f = self.anisotropic_form(0.05)
        self.assertTrue(is_zero(f.d2f_dgrad2[0, 1] - f.d2f_dgrad2[1, 0]))

    def test_convexity_threshold(self):
        # d2f/dp2 is positive definite iff W + W'' = W0 (1 - 15 e4 cos(4 theta)) > 0:
        # true for e4 = 0.05 < 1/15, false along theta = 0 for e4 = 0.1
        thetas = [k * math.pi / 32 for k in range(64)]
        for e4, convex in ((0.05, True), (0.1, False)):
            f = self.anisotropic_form(e4)
            smallest = min(
                min(float(ev) for ev in f.d2f_dgrad2.subs(self.values(e4, t)).eigenvals())
                for t in thetas)
            self.assertEqual(smallest > 0, convex, f"e4 = {e4}: smallest eigenvalue {smallest}")

    def test_finite_at_zero_gradient(self):
        f = self.anisotropic_form(0.05)
        values = self.values(0.05, 0.0, modulus=0.0)
        for expr in list(f.df_dgrad) + list(f.d2f_dgrad2):
            value = expr.subs(values)
            self.assertTrue(value.is_finite, f"{expr} is not finite at grad(phi) = 0")


class TestBuildErrors(unittest.TestCase):
    """
    Errors raised by build_variational_form on SymPy expressions.
    """

    def test_dimension_and_regularization(self):
        phi = sp.Symbol("phi")
        with self.assertRaises(VariationalError):
            build_variational_form(phi, ["phi"], 4)
        with self.assertRaises(VariationalError):
            build_variational_form(phi, ["phi"], 2, regularization=0.0)


class TestChecker(unittest.TestCase):
    """
    Consistency checks of the JSON attributes and expressions of variational coefficients.
    """

    def assertError(self, item, text):
        messages = errors(item)
        self.assertTrue(any(text in m for m in messages), f"'{text}' not in {messages}")

    def test_valid_coefficients(self):
        self.assertEqual(errors(coefficient("0.5*dot(grad(phi), grad(phi))")), [])
        self.assertEqual(errors(coefficient(
            f"0.5*(W0*{ANISOTROPY})**2*dot(grad(phi), grad(phi))",
            constants="(W0:1.0),(e4:0.05)", regularization=1.0e-8)), [])
        self.assertEqual(errors(coefficient("k*sdot(eta(1..2))/2 + c*eta1",
                                            variables="eta(1..2)", auxiliary_variables="c",
                                            constants="(k:1.0)")), [])

    def test_attributes(self):
        item = coefficient("dot(grad(phi), grad(phi))")
        del item["dimension"]
        self.assertError(item, "'dimension' (1, 2 or 3) must be defined")
        self.assertError(coefficient("dot(grad(phi), grad(phi))", dimension=4),
                         "'dimension' must be 1, 2 or 3")
        self.assertError(coefficient("dot(grad(phi), grad(phi))", dimension=True),
                         "'dimension' must be 1, 2 or 3")
        self.assertError(coefficient("dot(grad(phi), grad(phi))", regularization=-1.0),
                         "'regularization' must be a positive number")
        self.assertError(coefficient("dot(grad(phi), grad(phi))", variational="yes"),
                         "'variational' must be true or false")
        self.assertError(coefficient("dot(grad(phi), grad(phi))", gradient=True),
                         "cannot be both 'gradient' and 'variational'")

    def test_attributes_without_variational(self):
        item = {"class_name": "Test", "outputfile": "Test", "variables": "phi",
                "expression": "phi**2", "dimension": 2}
        self.assertError(item, "'dimension' is only used by a variational coefficient")

    def test_conditional_expression(self):
        item = coefficient("phi")
        del item["expression"]
        item["range_variable"] = "phi"
        item["expressions"] = [{"expression": "phi", "lower": 0.0, "upper": 1.0},
                               {"expression": "0.0", "default": True}]
        self.assertError(item, "not extended to variational coefficients")

    def test_grad_without_variational(self):
        item = {"class_name": "Test", "outputfile": "Test", "variables": "phi",
                "expression": "dot(grad(phi), grad(phi))"}
        self.assertError(item, "dot() and sdot() are only allowed")
        item["expression"] = "grad(phi)[0]"
        self.assertError(item, "grad() and n() are only allowed in a variational coefficient")

    def test_reserved_names(self):
        self.assertError(coefficient("dot(grad(phi), grad(phi))", auxiliary_variables="n"),
                         "is reserved in a variational coefficient")
        self.assertError(coefficient("dot(grad(phi), grad(phi))", constants="(grad_x:1.0)"),
                         "is reserved in a variational coefficient")

    def test_gradient_terms(self):
        self.assertError(coefficient("2*grad(phi)"), "must be used in dot(.,.) or with a component")
        self.assertError(coefficient("dot(phi, phi)"), "is neither grad(.) nor n(.)")
        self.assertError(coefficient("grad(phi)[2]"), "the component must be an integer between 0")
        self.assertError(coefficient("dot(grad(k), grad(phi))", constants="(k:1.0)"),
                         "is not a variable, an explicit or an auxiliary variable")
        self.assertError(coefficient("dot(grad(x), grad(phi))", coordinates="x,y"),
                         "is not a variable, an explicit or an auxiliary variable")
        self.assertError(coefficient("dot(grad(psi), grad(phi))"), "'psi' is used")

    def test_warning_without_gradient_of_a_variable(self):
        item = coefficient("phi**2*dot(grad(c), grad(c))", auxiliary_variables="c")
        self.assertEqual(errors(item), [])
        self.assertTrue(any("does not depend on the gradient" in m for m in warnings(item)))


class TestExistingFiles(unittest.TestCase):
    """
    Coefficient files of the SLOTH tests: no error with the checks of variational
    coefficients.
    """

    def test_no_error(self):
        files = []
        for path in sorted(TESTS_DIR.rglob("*.json")):
            try:
                data = json.loads(path.read_text())
            except (json.JSONDecodeError, UnicodeDecodeError):
                continue
            if isinstance(data, list) and data and all(
                    isinstance(item, dict) and "class_name" in item for item in data):
                files.append(path)
        self.assertTrue(files, f"no coefficient file found in {TESTS_DIR}")
        for path in files:
            with self.subTest(file=str(path.relative_to(TESTS_DIR))):
                found = [str(i) for i in check_file(path) if i.severity == "error"]
                self.assertEqual(found, [])


class TestCommandLine(unittest.TestCase):
    """
    GenerateCoefficient.py on a file with a variational coefficient: the check succeeds, the
    C++ generation stops with an error.
    """

    def run_script(self, *args, items=None):
        if items is None:
            items = [coefficient("0.5*dot(grad(phi), grad(phi))")]
        with tempfile.TemporaryDirectory() as tmp:
            json_file = Path(tmp) / "coefficients.json"
            json_file.write_text(json.dumps(items))
            env = dict(os.environ, PYTHONPATH=str(SCRIPTS_DIR), PYTHONDONTWRITEBYTECODE="1")
            result = subprocess.run([sys.executable, str(SCRIPTS_DIR / "GenerateCoefficient.py"),
                                     "-f", str(json_file), *args], cwd=tmp, env=env,
                                    capture_output=True, text=True)
            report = Path(tmp) / "coefficients_report.tex"
            result.report = report.read_text() if report.is_file() else None
            return result

    def test_check(self):
        result = self.run_script("--check")
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_report(self):
        item = coefficient(f"0.5*(W0*{ANISOTROPY})**2*dot(grad(phi), grad(phi))",
                           constants="(W0:1.0),(e4:0.05)")
        result = self.run_script("--report-only", items=[item])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIsNotNone(result.report)
        self.assertIn("variational (dimension 2)", result.report)
        self.assertIn("\\left|\\nabla \\phi\\right|", result.report)
        self.assertIn("\\mathbf{n}_{\\phi,x}", result.report)
        # derivatives written with n and |grad phi|, in the limit eps -> 0
        self.assertIn("\\frac{\\partial f}{\\partial (\\partial_{x} \\phi)}", result.report)
        self.assertIn("\\varepsilon \\to 0", result.report)
        self.assertNotIn("grad_phi", result.report)
        self.assertNotIn("Isotropic gradient term", result.report)

    def test_report_isotropic_and_data(self):
        items = [
            dict(coefficient("0.5*kappa*dot(grad(phi), grad(phi))", constants="(kappa:2.0)"),
                 class_name="Iso"),
            dict(coefficient("0.5*(k0 + k1*c)*dot(grad(phi), grad(phi))", variables="c",
                             auxiliary_variables="phi", constants="(k0:1.0),(k1:0.5)"),
                 class_name="Coupled"),
            dict(coefficient("0.5*k*sdot(eta(1..3))", variables="eta(1..3)", constants="(k:1.0)"),
                 class_name="Multi")]
        result = self.run_script("--report-only", items=items)
        self.assertEqual(result.returncode, 0, result.stderr)
        report = result.report
        self.assertIn("Isotropic gradient term", report)
        self.assertIn("\\lambda = \\kappa", report)
        self.assertIn("\\partial_{x} \\phi", report)
        self.assertIn("are data: they are not differentiated", report)
        self.assertIn("$9 \\times 9$ matrix (symmetric), not written", report)
        self.assertNotIn("varepsilon \\to 0", report)

    def test_generation_not_available(self):
        result = self.run_script()
        self.assertEqual(result.returncode, 1)
        self.assertIn("C++ generation not available for variational coefficients",
                      result.stderr)


if __name__ == "__main__":
    unittest.main()
