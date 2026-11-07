#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import subprocess
import json
from argparse import ArgumentParser, ArgumentTypeError, Namespace, RawTextHelpFormatter
from pathlib import Path
import os
import re
import sympy as sp


import shutil
import logging
import sys
from typing import List, Tuple, Optional

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s][%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)


# ======================
#    Constant Regex
# ======================

RANGE_PATTERN = re.compile(r'^([a-zA-Z_]\w*)\((\d+)\.\.(\d+)\)$')
SDOT_PATTERN = re.compile(r"(?:([-+]?\d*\.?\d+)\s*\*\s*)?sdot\((\w+)\((\d+)\.\.(\d+)\)\)")
DOT_PATTERN = re.compile(r"dot\(\s*(\w+)\s*,\s*(\w+)\s*\)")

def print_hessian(f, hessian, n):
    for i in range(n):
        for j in range(n):
            f.write(f"    hessian[{i*n + j}] = this->prefactor_ * ({sp.cxxcode(hessian[i,j])});\n")
            
def print_gradient(f, gradient, n):
    for k in range(n):
        f.write(f"    gradient[{k}] = this->prefactor_ * ({sp.cxxcode(gradient[k])});\n")

def get_coordinates(coordinate_names: str) -> List[str]:
    """
    Parse the optional 'coordinates' attribute, e.g. "X,Y" or "r,z".
    Names are mapped, in order, to the components of the physical point of the
    integration point (coordinates_[0], coordinates_[1], ...).
    """
    if not coordinate_names:
        return []
    names = [c.strip() for c in coordinate_names.split(',') if c.strip()]
    if len(names) > 3:
        raise ValueError("Error: at most three coordinates can be declared. Please check your data.")
    for c in names:
        if not re.match(r'^[a-zA-Z_]\w*$', c):
            raise ValueError(f"Error: invalid coordinate name '{c}'. Please check your data.")
    if len(set(names)) != len(names):
        raise ValueError("Error: duplicated coordinate names. Please check your data.")
    return names


def write_coordinates(f, coord_vars, class_name):
    """
    Write the local variables associated with the spatial coordinates.
    Coordinates are set by SLOTH (set_coordinates) before each evaluation.
    """
    if not coord_vars:
        return
    f.write(f"    if (this->coordinates_.size() < {len(coord_vars)}) {{\n")
    f.write(f"      throw std::runtime_error(\"{class_name}: spatial coordinates not set or "
            f"inconsistent with the dimension of the problem\");\n")
    f.write("    }\n")
    for i, c in enumerate(coord_vars):
        f.write(f"    double {c} = this->coordinates_[{i}];\n")


def check_bounds_consistency(expressions):
    for i in range(1, len(expressions)):
        _, _, _, _, upper_prev, upper_strict_prev = expressions[i-1]
        _, _, lower_curr, lower_strict_curr, _, _ = expressions[i]

        if (lower_curr is not None) and (upper_prev is not None) and (lower_curr < upper_prev):
            raise ValueError("Error: invalid bounds, intersection must not be empty. Please check your data.")

        if (lower_curr is not None) and (upper_strict_prev is not None) and (lower_curr < upper_strict_prev):
            raise ValueError("Error: invalid bounds, intersection must not be empty. Please check your data.")

        if (lower_strict_curr is not None) and (upper_prev is not None) and (lower_strict_curr < upper_prev):
            raise ValueError("Error: invalid bounds, intersection must not be empty. Please check your data.")

        if (lower_strict_curr is not None) and (upper_strict_prev is not None) and (lower_strict_curr < upper_strict_prev):
            raise ValueError("Error: invalid bounds, intersection must not be empty. Please check your data.")


def get_constants(constants):
    """
    Parse a string composed of couples (symbol, value), separated by semicolon
    Return a dict:
        - string of the constant
        - value of the constant
    """
    dict_constants = {t.replace("(","").replace(")","").split(':')[0]:t.replace("(","").replace(")","").split(':')[1] for t in constants.split(',')}
    # Translate R,NA,H,K constants as Physical constants in SLOTH
    for k,v in dict_constants.items():
        if v == 'K':
            dict_constants[k]="Physical::K"
        elif v == 'NA':
            dict_constants[k]="Physical::NA"
        elif v == 'H':
            dict_constants[k]="Physical::H"
        elif v == 'R':
            dict_constants[k]="Physical::R"
        elif v == 'T':
            dict_constants[k]="this->time_"
    return dict_constants


#  Case : a * sdot(x(1..j)) = a * (dot(x1,x1)+...+dot(xj,xj))
#  The expansion is always parenthesized so that any factor (numeric or symbolic, e.g.
#  kappa*sdot(x(1..j))) multiplies the whole sum.
def sdot_replace_fct(match):
    fact, var, start_str, end_str = match.groups()
    start, end = int(start_str), int(end_str)

    if start > end:
        raise ValueError("Error: invalid bounds in sdot expression. Please check your data.")

    expanded = "+".join(f"dot({var}{i},{var}{i})" for i in range(start, end+1))

    if fact:
        return f"{fact}*({expanded})"
    else:
        return f"({expanded})"
    
def expand_sdot(expr: str) -> str:
    """
    Parse a string composed of symbols and expand sdot .
    Return the expanded expression
    """
    return re.sub(SDOT_PATTERN, sdot_replace_fct, expr)


def print_header(f, list_expr, is_conditional_expression):
    f.write("/**\n")
    f.write(" *\n")
    f.write(" * @brief C++ function of the analytical expression\n")
    f.write(" *\n")

    # expression par défaut
    default_expr = list_expr[-1][0]
    if is_conditional_expression:
        f.write(" * Conditional definitions:\n")
        if_cond = "if "

        for i in range(len(list_expr) - 1):
            if i != 0:
                if_cond = "else if "
            expr, var, lower, lower_strict, upper, upper_strict = list_expr[i]

            cond_str = ""

            if lower is not None and upper is not None:
                cond_str = f"{lower} <= {var} <= {upper}"

            elif lower_strict is not None and upper_strict is not None:
                cond_str = f"{lower_strict} < {var} < {upper_strict}"

            elif lower is not None and upper_strict is not None:
                cond_str = f"{lower} <= {var} < {upper_strict}"

            elif lower_strict is not None and upper is not None:
                cond_str = f"{lower_strict} < {var} <= {upper}"

            f.write(f" *  {if_cond} ({cond_str})\n")
            f.write(f" *  {{ \n")
            f.write(f" *       F = {expr}\n")
            f.write(f" *  }}\n")
      
        f.write(f" *  else\n")
        f.write(f" *  {{ \n")
        f.write(f" *       F = {default_expr}\n")
        f.write(f" *  }}\n")
    else:
        f.write(f" *       F = {default_expr}\n")


    f.write(" */\n")
    
# -------------------
# -------------------

#  Case : a(i..j) in ai, ai+1, ..., aj.
def expand_ranges(expr: str,count: int = 0):
    """
    Parse a string composed of symbols and expand a(i..j) in ai, ai+1, ..., aj.
    Return a tuple:
        - string of expanded symbols
        - total number of variables expanded
    """
    tokens = [t.strip() for t in expr.split(',')]
    result = []

    for token in tokens:
        match = RANGE_PATTERN.match(token)
        if match:
            count += 1
            name, start, end = match.groups()
            start, end = int(start), int(end)

            if start > end:
                raise ValueError(f"Error: invalid bounds in the range {token}")

            for i in range(start, end + 1):
                result.append(f"{name}{i}")
        else:
            result.append(token)

    return ",".join(result), count

# -------------------
# -------------------
#  Case : dot(x,x) -> std::inner_product

def dot_to_inner_product(code: str) -> str:
    """
    Replace Sum(a[i]*b[i], (i, 0, n)) by std::inner_product(a.begin(), a.end(), b.begin(), 0)
    """

    def repl(match):
        var1 = match.group(1)
        var2 = match.group(2)
        return f"std::inner_product({var1}.begin(), {var1}.end(), {var2}.begin(), 0.0)"

    return DOT_PATTERN.sub(repl, code)


# -------------------
# -------------------

def sp_from_expr_with_sum(expression: str, local_dict, nb_var_expanded: int):
    """
    Convert string with Sum and IndexedBase into SymPy object.
    - Detect IndexedBase (variables)
    - Detect indexes of each summation
    - Fill locals_dict for sympify
    - Expand summation with .doit()
    - Delete []
    - Check if the number of expanded variables correspond to the number of IndexedBase
    """
    local_dict['Sum']= sp.Sum
    sum_vars = set(re.findall(r'([a-zA-Z_]\w*)\s*\[', expression))
    sum_indexes = set(re.findall(r'\[([a-zA-Z_]\w*)\]', expression))
    
    if len(sum_vars) != nb_var_expanded:
        raise Exception(f"The number of expanded variables must be equal to the number of IndexedBase. Please check json file")
    
    for v in sum_indexes:
        local_dict[v] = sp.symbols(v)
    for v in sum_vars:
        local_dict[v] = sp.IndexedBase(v)
    
    expr_tmp = sp.sympify(expression, locals=local_dict, rational=True)
    expr_tmp = sp.sympify(str(expr_tmp.doit()).replace('[', '').replace(']', ''), rational=True)
    
    return expr_tmp

# ==========================================================
#  Evaluation structure <class_name>Eval (host and device)
# ==========================================================
#
# In addition to the class <class_name>, a structure <class_name>Eval is written. It evaluates
# F, its gradient and its Hessian at one point with the same expressions, without any memory
# allocation, and can be called both on the host and on the device (MFEM_HOST_DEVICE), for
# instance inside mfem::forall (see Coefficients/EvalAtPoints.hpp).
#
# Names of the arguments of <class_name>Eval::operator(): an expression using one of them as a
# variable, a constant or a coordinate is not written as an evaluation structure.
EVAL_RESERVED_NAMES = {"values", "explicit_values", "auxiliary_values", "coordinates", "time",
                       "F", "gradient", "hessian", "prefactor_"}


def get_bounds(bounds, variable_names, constant_names, class_name):
    """
    Parse the optional 'bounds' attribute: bounds of variables, e.g.
        "bounds": {"x": ["eps", "1 - eps"], "eta(1..2)": [0.0, null]}
    A key is a variable (unknown, explicit or auxiliary variable) or a range of variables
    eta(1..2). A bound is a number, an expression of the constants of the coefficient, or null
    (no bound on this side).

    The value of a bounded variable is clamped to its bounds before the evaluation of F, of its
    gradient and of its Hessian, which are thus evaluated at the clamped value (outside the
    bounds, the derivatives are those at the bound, not the zero derivatives of the clamped
    function).

    Return a dict: variable name -> (lower bound, upper bound), each one a C++ expression or None.
    """
    if not bounds:
        return {}
    if not isinstance(bounds, dict):
        raise ValueError(f"Error ({class_name}): 'bounds' must be a dictionary "
                         "{variable: [lower, upper]}. Please check your data.")

    def to_cpp(bound):
        if bound is None:
            return None
        if isinstance(bound, bool) or not isinstance(bound, (int, float, str)):
            raise ValueError(f"Error ({class_name}): a bound must be a number, an expression of "
                             f"the constants or null, not {bound!r}. Please check your data.")
        if isinstance(bound, (int, float)):
            return repr(float(bound))
        if "**" in bound or "^" in bound:
            raise ValueError(f"Error ({class_name}): powers are not allowed in the bound "
                             f"'{bound}'. Please check your data.")
        try:
            symbols = {str(x) for x in sp.sympify(bound).free_symbols}
        except (sp.SympifyError, SyntaxError, TypeError) as e:
            raise ValueError(f"Error ({class_name}): invalid bound '{bound}'. "
                             "Please check your data.") from e
        unknown = symbols - set(constant_names)
        if unknown:
            raise ValueError(f"Error ({class_name}): the bound '{bound}' uses {sorted(unknown)}, "
                             "which are not constants of the coefficient. Please check your data.")
        bound = bound.strip()
        if re.fullmatch(r"[A-Za-z_]\w*", bound):
            return bound
        return f"({bound})"

    result = {}
    for key, value in bounds.items():
        match = RANGE_PATTERN.match(key.strip())
        if match:
            name, start, end = match.groups()
            names = [f"{name}{i}" for i in range(int(start), int(end) + 1)]
        else:
            names = [key.strip()]
        if not isinstance(value, (list, tuple)) or len(value) != 2:
            raise ValueError(f"Error ({class_name}): the bounds of '{key}' must be a list "
                             "[lower, upper]. Please check your data.")
        lower, upper = value
        if lower is None and upper is None:
            raise ValueError(f"Error ({class_name}): at least one bound of '{key}' must be "
                             "defined. Please check your data.")
        if isinstance(lower, (int, float)) and isinstance(upper, (int, float)) and lower > upper:
            raise ValueError(f"Error ({class_name}): lower bound greater than upper bound for "
                             f"'{key}'. Please check your data.")
        for v in names:
            if v not in variable_names:
                raise ValueError(f"Error ({class_name}): '{v}' in 'bounds' is not a variable of "
                                 "the coefficient. Please check your data.")
            if v in result:
                raise ValueError(f"Error ({class_name}): bounds of '{v}' defined twice. "
                                 "Please check your data.")
            result[v] = (to_cpp(lower), to_cpp(upper))
    return result


def bounded_value(value, variable_bounds):
    """
    C++ expression of a variable read from value and clamped to its bounds (if any).
    """
    if variable_bounds is None:
        return value
    lower, upper = variable_bounds
    if lower is not None and upper is not None:
        return f"fmin({upper}, fmax({lower}, {value}))"
    if lower is not None:
        return f"fmax({lower}, {value})"
    return f"fmin({upper}, {value})"


def condition_to_cpp(item):
    """
    C++ condition of one piece of a conditional expression.
    item = (expression, variable, lower, lower_strict, upper, upper_strict)
    """
    _, var, lower, lower_strict, upper, upper_strict = item
    left = f"({lower} <= {var})" if lower is not None else f"({lower_strict} < {var})"
    right = f"({var} <= {upper})" if upper is not None else f"({var} < {upper_strict})"
    return f"{left} && {right}"


def write_eval_expressions(f, expr, gradient, hessian, n, indent):
    """
    Write F, its gradient and its Hessian (each one only if requested by a non-null pointer).
    The Hessian is symmetric: its lower part is copied from its upper part.
    """
    f.write(f"{indent}if (F != nullptr) {{\n")
    f.write(f"{indent}  *F = prefactor_ * ({sp.ccode(expr)});\n")
    f.write(f"{indent}}}\n")
    f.write(f"{indent}if (gradient != nullptr) {{\n")
    for k in range(n):
        f.write(f"{indent}  gradient[{k}] = prefactor_ * ({sp.ccode(gradient[k])});\n")
    f.write(f"{indent}}}\n")
    f.write(f"{indent}if (hessian != nullptr) {{\n")
    for i in range(n):
        for j in range(n):
            if j < i:
                f.write(f"{indent}  hessian[{i*n + j}] = hessian[{j*n + i}];\n")
            else:
                f.write(f"{indent}  hessian[{i*n + j}] = prefactor_ * ({sp.ccode(hessian[i, j])});\n")
    f.write(f"{indent}}}\n")


def write_eval_struct(cpp_file, class_name, list_expr, list_gradient, list_hessian, spvars,
                      explicit_vars, auxiliary_vars, dict_constants, coord_vars,
                      is_conditional_expression, bounds=None):
    """
    Write the structure <class_name>Eval in cpp_file (see the comment above EVAL_RESERVED_NAMES).
    bounds: bounds of variables (see get_bounds), applied when the variables are read.

    Hessian stored in a vector: H(i,j) -> H(i*n+j), n number of variables.
    Constants are members of the structure (copied with it on the device), except the time (T),
    given as an argument.
    Return True if the structure is written.
    """
    n = len(spvars)
    names = set(str(v) for v in spvars) | set(str(v) for v in explicit_vars) \
        | set(str(v) for v in auxiliary_vars) | set(dict_constants.keys()) | set(coord_vars)
    clash = names & EVAL_RESERVED_NAMES
    bounds = bounds or {}
    if clash and bounds:
        raise ValueError(f"Error ({class_name}): names {sorted(clash)} are reserved for the "
                         "evaluation structure, required by 'bounds'. Please rename them.")
    if clash:
        logging.warning(f"{class_name}Eval not generated: names {sorted(clash)} are reserved "
                        f"for the evaluation structure. The class {class_name} holds its own "
                        "expressions.")
        return False

    time_constants = [k for k, v in dict_constants.items() if v == "this->time_"]
    member_constants = {k: v for k, v in dict_constants.items() if v != "this->time_"}

    with open(cpp_file, "a") as f:
        f.write("/**\n")
        f.write(f" * @brief Evaluation of {class_name} at one point, on the host or on the device\n")
        f.write(" *\n")
        f.write(f" * Same expressions as {class_name} (F, gradient and Hessian), without memory\n")
        f.write(" * allocation. Each output is computed only if its pointer is not null.\n")
        if bounds:
            f.write(" *\n")
            f.write(" * Bounded variables (F, gradient and Hessian evaluated at the clamped values):\n")
            for v, (lower, upper) in bounds.items():
                lo = lower if lower is not None else "-inf"
                up = upper if upper is not None else "+inf"
                f.write(f" *   {v} in [{lo}, {up}]\n")
        f.write(" * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)\n")
        f.write(" */\n")
        f.write(f"struct {class_name}Eval {{\n")
        f.write(f"  static constexpr int nb_variables = {n};\n")
        f.write(f"  static constexpr int nb_explicit_variables = {len(explicit_vars)};\n")
        f.write(f"  static constexpr int nb_auxiliary_variables = {len(auxiliary_vars)};\n")
        f.write(f"  static constexpr int nb_coordinates = {len(coord_vars)};\n\n")
        f.write("  double prefactor_ = 1.0;\n")
        for k, v in member_constants.items():
            f.write(f"  double {k} = {v};\n")
        f.write("\n")
        f.write("  MFEM_HOST_DEVICE void operator()([[maybe_unused]] const double* values,\n")
        f.write("                                   [[maybe_unused]] const double* explicit_values,\n")
        f.write("                                   [[maybe_unused]] const double* auxiliary_values,\n")
        f.write("                                   [[maybe_unused]] const double* coordinates,\n")
        f.write("                                   [[maybe_unused]] const double time, double* F,\n")
        f.write("                                   double* gradient, double* hessian) const {\n")
        if bounds:
            f.write("    // Bounded variables: clamped to their bounds\n")
        for i, v in enumerate(spvars):
            value = bounded_value(f"values[{i}]", bounds.get(str(v)))
            f.write(f"    [[maybe_unused]] const double {v} = {value};\n")
        for i, v in enumerate(explicit_vars):
            value = bounded_value(f"explicit_values[{i}]", bounds.get(str(v)))
            f.write(f"    [[maybe_unused]] const double {v} = {value};\n")
        for i, v in enumerate(auxiliary_vars):
            value = bounded_value(f"auxiliary_values[{i}]", bounds.get(str(v)))
            f.write(f"    [[maybe_unused]] const double {v} = {value};\n")
        for i, c in enumerate(coord_vars):
            f.write(f"    [[maybe_unused]] const double {c} = coordinates[{i}];\n")
        for k in time_constants:
            f.write(f"    [[maybe_unused]] const double {k} = time;\n")

        if is_conditional_expression:
            for i in range(len(list_expr) - 1):
                keyword = "if" if i == 0 else "} else if"
                f.write(f"    {keyword} ({condition_to_cpp(list_expr[i])}) {{\n")
                write_eval_expressions(f, list_expr[i][0], list_gradient[i][0],
                                       list_hessian[i][0], n, "      ")
            f.write("    } else {\n")
            write_eval_expressions(f, list_expr[-1][0], list_gradient[-1][0],
                                   list_hessian[-1][0], n, "      ")
            f.write("    }\n")
        else:
            write_eval_expressions(f, list_expr[-1][0], list_gradient[-1][0],
                                   list_hessian[-1][0], n, "    ")
        f.write("  }\n")
        f.write("};\n\n")
    return True


#######################################################################
#######################################################################


def prepare_output_file(output_file, is_gradient_coefficient):

    filename = Path(output_file).resolve()

    with open(output_file, "w") as f:
        # Header 
        f.write("""/**
 *
 * Copyright CEA (C) 2026
 *
 * This file is part of SLOTH.
 *
 * SLOTH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLOTH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */\n
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>\n

#include "Options/PhysicalPropertiesOptions.hpp"\n  
#include "Coefficients/FunctionCoefficient.hpp"\n  
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]\n

#pragma once\n
""")
        
    
def write_coordinates_check(f, coord_vars, class_name):
    """
    Write the check of the spatial coordinates set by SLOTH (set_coordinates) before each
    evaluation.
    """
    if not coord_vars:
        return
    f.write(f"    if (this->coordinates_.size() < {len(coord_vars)}) {{\n")
    f.write(f"      throw std::runtime_error(\"{class_name}: spatial coordinates not set or "
            f"inconsistent with the dimension of the problem\");\n")
    f.write("    }\n")


def write_class_using_eval(cpp_file, class_name, list_expr, n, coord_vars,
                           is_conditional_expression):
    """
    Write the class <class_name>: F(), GradientF() and HessianF() call the evaluation structure
    <class_name>Eval (written before the class), which holds the expressions.
    """
    args = ("const std::span<const double>&, const std::span<const double>&, "
            "const std::span<const double>&, const unsigned int dimension")
    lambda_head = ("  auto func = [&](const std::span<const double>& input_vector,\n"
                   "                  const std::span<const double>& exp_input_vector,\n"
                   "                  const std::span<const double>& auxiliary_vector,\n"
                   "                  [[maybe_unused]] const unsigned int dimension) {\n")

    def eval_call(outputs):
        return ("    this->eval_(input_vector.data(), exp_input_vector.data(), "
                "auxiliary_vector.data(),\n"
                f"                this->coordinates_.data(), this->time_, {outputs});\n")

    with open(cpp_file, "a") as f:
        print_header(f, list_expr, is_conditional_expression)

        # Class
        f.write(f"class {class_name} : public FunctionCoefficient {{\n")
        f.write(" private:\n")
        f.write(f"  {class_name}Eval eval_;  ///< F, gradient and Hessian (host and device)\n\n")
        f.write(" protected:\n")
        f.write(f"  std::function<double({args})> F() final;\n")
        f.write(f"  std::function<std::vector<double>({args})> GradientF() final;\n")
        f.write(f"  std::function<std::vector<double>({args})> HessianF() final;\n\n")
        f.write(" public:\n")
        f.write(f"  {class_name}() = default;\n")
        f.write(f"  explicit {class_name}(const double prefactor) {{ eval_.prefactor_ = prefactor; }}\n")
        f.write(f"  virtual ~{class_name}() = default;\n")
        f.write(f"  /// Evaluation structure (constants, prefactor), e.g. for eval_at_points\n")
        f.write(f"  const {class_name}Eval& get_eval() const {{ return eval_; }}\n")
        if coord_vars:
            f.write("  bool uses_coordinates() const final { return true; }\n")
        f.write("};\n\n")

        # F
        f.write("/**\n")
        f.write(f" * @brief Value of {class_name} (see {class_name}Eval)\n")
        f.write(" */\n")
        f.write(f"std::function<double({args})> {class_name}::F() {{\n")
        f.write(lambda_head)
        write_coordinates_check(f, coord_vars, class_name)
        f.write("    double value = 0.0;\n")
        f.write(eval_call("&value, nullptr, nullptr"))
        f.write("    return value;\n")
        f.write("  };\n")
        f.write("  return func;\n")
        f.write("}\n\n")

        # Gradient
        f.write("/**\n")
        f.write(f" * @brief Gradient of {class_name} (see {class_name}Eval)\n")
        f.write(" */\n")
        f.write(f"std::function<std::vector<double>({args})> {class_name}::GradientF() {{\n")
        f.write(lambda_head)
        write_coordinates_check(f, coord_vars, class_name)
        f.write(f"    std::vector<double> gradient({n});\n")
        f.write(eval_call("nullptr, gradient.data(), nullptr"))
        f.write("    return gradient;\n")
        f.write("  };\n")
        f.write("  return func;\n")
        f.write("}\n\n")

        # Hessian
        f.write("/**\n")
        f.write(f" * @brief Hessian of {class_name} (see {class_name}Eval)\n")
        f.write(" * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)\n")
        f.write(" */\n")
        f.write(f"std::function<std::vector<double>({args})> {class_name}::HessianF() {{\n")
        f.write(lambda_head)
        write_coordinates_check(f, coord_vars, class_name)
        f.write(f"    std::vector<double> hessian({n * n});\n")
        f.write(eval_call("nullptr, nullptr, hessian.data()"))
        f.write("    return hessian;\n")
        f.write("  };\n")
        f.write("  return func;\n")
        f.write("}\n\n")


def write_class_with_expressions(cpp_file, class_name, list_expr, list_gradient, list_hessian,
                                 spvars, n, explicit_vars, auxiliary_vars, dict_constants,
                                 coord_vars, has_explicit_variables, has_auxiliary_variables,
                                 has_constants, has_coordinates, is_gradient_coefficient,
                                 is_conditional_expression):
    """
    Write the class <class_name> holding its own expressions in F(), GradientF() and HessianF().
    Used when the evaluation structure <class_name>Eval is not generated (gradient expressions,
    reserved names).
    """
    with open(cpp_file, "a") as f:

        print_header(f, list_expr, is_conditional_expression)
 
        # Class
        f.write(f"class {class_name} : public FunctionCoefficient {{\n")
        f.write(" private:\n")
        f.write("  double prefactor_;\n")
        f.write(" protected:\n")
        f.write("  std::function<double(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> F() final;\n")
        f.write("  std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> GradientF() final;\n")
        f.write("  std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&,const std::span<const double>&, const unsigned int dimension)> HessianF() final;\n\n")
        f.write(" public:\n")
        f.write(f"  {class_name}() : prefactor_(1.0) {{}}\n")
        f.write(f" explicit {class_name}(const double prefactor): prefactor_(prefactor) {{}}\n")
        f.write(f"  virtual ~{class_name}() = default;\n")
        if has_coordinates:
            f.write("  bool uses_coordinates() const final { return true; }\n")
        f.write("};\n\n")


        f.write("/**\n")
        f.write(" *\n")
        f.write(f" * @brief C++ function of the expression\n")
        f.write(" *\n")

        # Expression par défaut
        default_expr = list_expr[-1][0]
        f.write(f" *\n")
        if is_conditional_expression:
            if_cond = "if "
            for i in range(0, len(list_expr)-1):
                if i != 0:
                    if_cond = "else if "
                expr, var, lower, lower_strict, upper, upper_strict = list_expr[i]

                cond_str = ""

                # <= && >=
                if lower is not None and upper is not None:
                    cond_str = f"{lower} <= {var} <= {upper}"

                # < && >
                elif lower_strict is not None and upper_strict is not None:
                    cond_str = f"{lower_strict} < {var} < {upper_strict}"

                # <= && >
                elif lower is not None and upper_strict is not None:
                    cond_str = f"{lower} <= {var} < {upper_strict}"

                # < && >=
                elif lower_strict is not None and upper is not None:
                    cond_str = f"{lower_strict} < {var} <= {upper}"

                f.write(f" *  {if_cond} ({cond_str})\n")
                f.write(f" *  {{\n")
                f.write(f" *       F = {sp.cxxcode(expr)}\n")
                f.write(f" *  }}\n")
            f.write(f" *  else\n")
            f.write(f" *  {{\n")
            f.write(f" *       F = {sp.cxxcode(default_expr)}\n")
            f.write(f" *  }}\n")  
        
        f.write(f" * @return std::function<double(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)>\n")
        f.write(f" */\n ")  
        
        # Fonction F()
        f.write(f"std::function<double(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> {class_name}::F() {{\n")

        if(not is_gradient_coefficient):        
            if(has_auxiliary_variables):      
                if(has_explicit_variables):
                    f.write("  auto func = [&](const std::span<const double>& input_vector, const std::span<const double>& exp_input_vector, const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {\n")
                else:
                    f.write("  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&,const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {\n")
            else:     
                if(has_explicit_variables):
                    f.write("  auto func = [&](const std::span<const double>& input_vector, const std::span<const double>& exp_input_vector,[[maybe_unused]] const std::span<const double>&,  [[maybe_unused]] const unsigned int dimension) {\n")
                else:
                    f.write("  auto func = [&](const std::span<const double>& input_vector,[[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const unsigned int dimension) {\n")
            
            for i, v in enumerate(spvars):
                f.write(f"    double {v} = input_vector[{i}];\n")
            if(has_auxiliary_variables):
                for i, v in enumerate(auxiliary_vars):
                    f.write(f"    double {v} = auxiliary_vector[{i}];\n")
            if(has_explicit_variables):
                for i, v in enumerate(explicit_vars):
                    f.write(f"    double {v} = exp_input_vector[{i}];\n")
            
            if(has_constants):
                for const in dict_constants:
                    f.write(f"    double {const} = {dict_constants[const]};\n")
            write_coordinates(f, coord_vars, class_name)
   
            #====================================================
            #====================================================
            # single expression or the default one in case "conditional expression"
            expr_0 = list_expr[-1][0]
            if is_conditional_expression:
                f.write(f"    double F = 0.0;\n")
                if_cond = "if "
                for i in range(0,len(list_expr)-1):
                    if i != 0:
                        if_cond = "else if "
                    expr = list_expr[i][0]
                    # ===============
                    # <= && >=
                    if list_expr[i][2] is not None and list_expr[i][4] is not None:
                        f.write(f"  {if_cond}(({list_expr[i][2]}<={list_expr[i][1]}) && ({list_expr[i][1]}<={list_expr[i][4]}))  \n")
                        f.write("{ \n")
                        f.write(f"    F = {sp.cxxcode(expr)};\n")
                        f.write("} \n")
                    # ===============
                    # < && >
                    if list_expr[i][3] is not None and list_expr[i][5] is not None:
                        f.write(f"   {if_cond}(({list_expr[i][3]}<{list_expr[i][1]}) && ({list_expr[i][1]}<{list_expr[i][5]}))  \n")
                        f.write("{ \n")
                        f.write(f"    F = {sp.cxxcode(expr)};\n")
                        f.write("} \n")
                    # ===============
                    # <= && >
                    if list_expr[i][2] is not None and list_expr[i][5] is not None:
                        f.write(f"   {if_cond}(({list_expr[i][2]}<={list_expr[i][1]}) && ({list_expr[i][1]}<{list_expr[i][5]}))  \n")
                        f.write("{ \n")
                        f.write(f"    F = {sp.cxxcode(expr)};\n")
                        f.write("} \n")
                    # ===============
                    # < && >=
                    if list_expr[i][3] is not None and list_expr[i][4] is not None:
                        f.write(f"   {if_cond}(({list_expr[i][3]}<{list_expr[i][1]}) && ({list_expr[i][1]}<={list_expr[i][4]}))  \n")
                        f.write("{ \n")
                        f.write(f"    F = {sp.cxxcode(expr)};\n")
                        f.write("} \n")
                f.write(f" else \n")
                f.write("{ \n")
                f.write(f"    F = {sp.cxxcode(expr_0)};\n")
                f.write("} \n")
            else:
                f.write(f" double   F = {sp.cxxcode(expr_0)};\n")

                
        else:        
            if(has_auxiliary_variables):      
                if(has_explicit_variables):
                    f.write("  auto func = [&](const std::span<const double>& input_vector, const std::span<const double>& exp_input_vector, const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {\n")
                else:
                    f.write("  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&,const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {\n")
            else:     
                if(has_explicit_variables):
                    f.write("  auto func = [&](const std::span<const double>& input_vector, const std::span<const double>& exp_input_vector,[[maybe_unused]] const std::span<const double>&,  [[maybe_unused]] const unsigned int dimension) {\n")
                else:
                    f.write("  auto func = [&](const std::span<const double>& input_vector,[[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const unsigned int dimension) {\n")
            
            for i, v in enumerate(spvars):
                f.write(f"    std::vector<double> {v};\n")
                f.write(f"    for(unsigned int i=0;i<dimension;i++) {v}.push_back(input_vector[{i}*dimension+i]);\n")

            if(has_auxiliary_variables):
                for i, v in enumerate(auxiliary_vars):
                    f.write(f"    std::vector<double> {v};\n")
                    f.write(f"    for(unsigned int i=0;i<dimension;i++) {v}.push_back(auxiliary_vector[{i}*dimension+i]);\n")

            if(has_explicit_variables):
                for i, v in enumerate(explicit_vars):
                    f.write(f"    std::vector<double> {v};\n")
                    f.write(f"    for(unsigned int i=0;i<dimension;i++) {v}.push_back(exp_input_vector[{i}*dimension+i]);\n")

            if(has_constants):
                for const in dict_constants:
                    f.write(f"    double {const} = {dict_constants[const]};\n")
            expr_0 = list_expr[-1][0]
            f.write(f"    double F = {dot_to_inner_product(str(sp.N(expr_0)))};\n")
            if is_conditional_expression:
                raise ValueError("Error: conditional expression not extended to gradient expression.")  

   
        f.write("    return this->prefactor_ * F;\n")
        f.write("  };\n")
        f.write("  return func;\n")
        f.write("}\n\n")

        f.write(f"""/**
 *
 * @brief Gradient
 * 
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const unsigned int dimension)> 
 */
""")
        # Gradient
        f.write(f"std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> {class_name}::GradientF() {{\n")


        if(not is_gradient_coefficient):   
            if(has_auxiliary_variables):      
                if(has_explicit_variables):
                    f.write("  auto func = [&](const std::span<const double>& input_vector, const std::span<const double>& exp_input_vector, const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {\n")
                else:
                    f.write("  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&,const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {\n")
            else:     
                if(has_explicit_variables):
                    f.write("  auto func = [&](const std::span<const double>& input_vector, const std::span<const double>& exp_input_vector,[[maybe_unused]] const std::span<const double>&,  [[maybe_unused]] const unsigned int dimension) {\n")
                else:
                    f.write("  auto func = [&](const std::span<const double>& input_vector,[[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const unsigned int dimension) {\n")
            
            for i, v in enumerate(spvars):
                f.write(f"    double {v} = input_vector[{i}];\n")
                
            if(has_auxiliary_variables):
                for i, v in enumerate(auxiliary_vars):
                    f.write(f"    double {v} = auxiliary_vector[{i}];\n")
                    
            if(has_explicit_variables):
                for i, v in enumerate(explicit_vars):
                    f.write(f"    double {v} = exp_input_vector[{i}];\n")
                    
            if(has_constants):
                for const in dict_constants:
                    f.write(f"    double {const} = {dict_constants[const]};\n")
            write_coordinates(f, coord_vars, class_name)
                    
            f.write(f"    std::vector<double> gradient({n});\n")

            #====================================================
            #====================================================
            # single gradient or the default one in case "conditional expression"
            gradient_0 = list_gradient[-1][0]

            if is_conditional_expression:
                if_cond = "if "
                for i in range(0,len(list_gradient)-1):
                    if i != 0:
                        if_cond = "else if "
                    gradient = list_gradient[i][0]
                    # ===============
                    # <= && >=
                    if list_gradient[i][2] is not None and list_gradient[i][4] is not None:
                        f.write(f"  {if_cond}(({list_gradient[i][2]}<={list_gradient[i][1]}) && ({list_gradient[i][1]}<={list_gradient[i][4]}))  \n")
                        f.write("{ \n")
                        print_gradient(f, gradient, n)
                        f.write("} \n")
                    # ===============
                    # < && >
                    if list_gradient[i][3] is not None and list_gradient[i][5] is not None:
                        f.write(f"  {if_cond}(({list_gradient[i][3]}<{list_gradient[i][1]}) && ({list_gradient[i][1]}<{list_gradient[i][5]}))  \n")
                        f.write("{ \n")
                        print_gradient(f, gradient, n)
                        f.write("} \n")
                    # ===============
                    # <= && >
                    if list_gradient[i][2] is not None and list_gradient[i][5] is not None:
                        f.write(f"  {if_cond}(({list_gradient[i][2]}<={list_gradient[i][1]}) && ({list_gradient[i][1]}<{list_gradient[i][5]}))  \n")
                        f.write("{ \n")
                        print_gradient(f, gradient, n)
                        f.write("} \n")
                    # ===============
                    # < && >=
                    if list_gradient[i][3] is not None and list_gradient[i][4] is not None:
                        f.write(f"  {if_cond}(({list_gradient[i][3]}<{list_gradient[i][1]}) && ({list_gradient[i][1]}<={list_gradient[i][4]}))  \n")
                        f.write("{ \n")
                        print_gradient(f, gradient, n)
                        f.write("} \n")
                f.write(f" else \n")
                f.write("{ \n")
                print_gradient(f, gradient_0, n)
                f.write("} \n")
            else:
                print_gradient(f, gradient_0, n)

        else: 
            f.write("  auto func = [&]([[maybe_unused]] const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const unsigned int dimension) {\n")
            f.write(f"    std::vector<double> gradient({n},0.0);\n")

        f.write("    return gradient;\n")
        f.write("  };\n")
        f.write("  return func;\n")
        f.write("}\n\n")


        f.write(f"""/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 * 
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> 
 */
""")
        # HessianF()
        f.write(f"std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> {class_name}::HessianF() {{\n")
        if(not is_gradient_coefficient): 
            if(has_auxiliary_variables):      
                if(has_explicit_variables):
                    f.write("  auto func = [&](const std::span<const double>& input_vector, const std::span<const double>& exp_input_vector, const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {\n")
                else:
                    f.write("  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&,const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {\n")
            else:     
                if(has_explicit_variables):
                    f.write("  auto func = [&](const std::span<const double>& input_vector, const std::span<const double>& exp_input_vector,[[maybe_unused]] const std::span<const double>&,  [[maybe_unused]] const unsigned int dimension) {\n")
                else:
                    f.write("  auto func = [&](const std::span<const double>& input_vector,[[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const unsigned int dimension) {\n")
            
            for i, v in enumerate(spvars):
                f.write(f"    double {v} = input_vector[{i}];\n")
                
            if(has_auxiliary_variables):
                for i, v in enumerate(auxiliary_vars):
                    f.write(f"    double {v} = auxiliary_vector[{i}];\n")
                    
            if(has_explicit_variables):
                for i, v in enumerate(explicit_vars):
                    f.write(f"    double {v} = exp_input_vector[{i}];\n")
                
            if(has_constants):
                for const in dict_constants:
                    f.write(f"    double {const} = {dict_constants[const]};\n")
            write_coordinates(f, coord_vars, class_name)
                    
            f.write(f"    std::vector<double> hessian({n*n});\n")
       
            #====================================================
            #====================================================
            # single hessian or the default one in case "conditional expression"
            hessian_0 = list_hessian[-1][0]
                        
            if is_conditional_expression:
                if_cond = "if "
                for i in range(0,len(list_hessian)-1):
                    if i != 0:
                        if_cond = "else if "
                    hessian = list_hessian[i][0]
                    # ===============
                    # <= && >=
                    if list_hessian[i][2] is not None and list_hessian[i][4] is not None:
                        f.write(f"  {if_cond}(({list_hessian[i][2]}<={list_hessian[i][1]}) && ({list_hessian[i][1]}<={list_hessian[i][4]}))  \n")
                        f.write("{ \n")
                        print_hessian(f, hessian, n)
                        f.write("} \n")
                    # ===============
                    # < && >
                    if list_hessian[i][3] is not None and list_hessian[i][5] is not None:
                        f.write(f"  {if_cond}(({list_hessian[i][3]}<{list_hessian[i][1]}) && ({list_hessian[i][1]}<{list_hessian[i][5]}))  \n")
                        f.write("{ \n")
                        print_hessian(f, hessian, n)
                        f.write("} \n")
                    # ===============
                    # <= && >
                    if list_hessian[i][2] is not None and list_hessian[i][5] is not None:
                        f.write(f" {if_cond}(({list_hessian[i][2]}<={list_hessian[i][1]}) && ({list_hessian[i][1]}<{list_hessian[i][5]}))  \n")
                        f.write("{ \n")
                        print_hessian(f, hessian, n)
                        f.write("} \n")
                    # ===============
                    # < && >=
                    if list_hessian[i][3] is not None and list_hessian[i][4] is not None:
                        f.write(f"  {if_cond}(({list_hessian[i][3]}<{list_hessian[i][1]}) && ({list_hessian[i][1]}<={list_hessian[i][4]}))  \n")
                        f.write("{ \n")
                        print_hessian(f, hessian, n)
                        f.write("} \n")
                f.write(f" else \n")
                f.write("{ \n")
                print_hessian(f, hessian_0, n)
                f.write("} \n")
            else:
                print_hessian(f, hessian_0, n)


        else:
            f.write("  auto func = [&]([[maybe_unused]] const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&,  [[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const unsigned int dimension) {\n")
            f.write(f"    std::vector<double> hessian({n*n},0.0);\n")

        f.write("    return hessian;\n")
        f.write("  };\n")
        f.write("  return func;\n")
        f.write("}\n")


def generate_class_with_functions(list_expr_str1, var_names, explicit_var_names, auxiliary_var_names, constants, class_name, output_file, is_gradient_coefficient, is_conditional_expression, coordinate_names="", bounds=None):
    
    # Get constants
    has_constants = False
    if(len(constants)>0):
        dict_constants = get_constants(constants)
        has_constants = True
        constants_names=','.join(list(dict_constants.keys()))+ ","
        constants_vars = sp.symbols(constants_names)
    
    
    # Expand variables with a range if needed
    nb_var_expanded = 0
    var_names, nb_var_expanded = expand_ranges(var_names, nb_var_expanded)
    var_names+=","
    spvars = sp.symbols(var_names)
    n = len(spvars)
    has_auxiliary_variables = False
    if(len(auxiliary_var_names)>0):
        auxiliary_var_names, nb_var_expanded = expand_ranges(auxiliary_var_names, nb_var_expanded)
        has_auxiliary_variables = True
        auxiliary_var_names+=","
        auxiliary_vars = sp.symbols(auxiliary_var_names)
        
    has_explicit_variables = False
    if(len(explicit_var_names)>0):
        explicit_var_names, nb_var_expanded = expand_ranges(explicit_var_names, nb_var_expanded)
        has_explicit_variables = True
        explicit_var_names+=","
        explicit_vars = sp.symbols(explicit_var_names)

    # Spatial coordinates (optional): not differentiated (they are not unknowns)
    coord_vars = get_coordinates(coordinate_names)
    has_coordinates = len(coord_vars) > 0
    if has_coordinates:
        if is_gradient_coefficient:
            raise ValueError("Error: coordinates are not extended to gradient expression.")
        other_names = set(str(v) for v in spvars)
        if has_auxiliary_variables:
            other_names |= set(str(v) for v in auxiliary_vars)
        if has_explicit_variables:
            other_names |= set(str(v) for v in explicit_vars)
        if has_constants:
            other_names |= set(dict_constants.keys())
        clash = other_names & set(coord_vars)
        if clash:
            raise ValueError(f"Error: coordinate names {sorted(clash)} are already used as variables "
                             "or constants. Please check your data.")

    # translate sdot contributions before managing expressions with sympy    
    # pattern = r"(?:([-+]?\d*\.?\d+)\s*\*\s*)?sdot\((\w+)\((\d+)\.\.(\d+)\)\)"

    # ==========================================================
    # ========================================================== 
    # List of expressions and associated derivatives 
    list_expr = []
    list_gradient = []
    list_hessian = []
    # Get all expressions even for those of type gradient (simpler approach)
    for expr_str1 in list_expr_str1:
        
        expr_str = expand_sdot(expr_str1[0])
        locals_dict={}
        
        # Check if the expression is of type gradient 
        has_dot = bool(re.search(r"\bdot\s*\(", expr_str))
        if has_dot and (not is_gradient_coefficient):
            raise ValueError(f"Analytical expression contains at least a dot(...) term but "
                            "is not declared as a gradient expression."
                            "Expression of type gradient are differentiated.")

        if(is_gradient_coefficient):
            dot = sp.Function('dot')
            locals_dict["dot"]=dot
            
        # Check presence of Sum terms in the expression
        bad_sum_term = bool(re.search(r"\bsum\s*\(", expr_str))
        if bad_sum_term:
            raise ValueError(f"Analytical expression contains at least a sum term but is incorrectly written."
                            "Summation is defined Sum(...)")
            
        has_sum = "Sum(" in expr_str
        if has_sum:
            expr_tmp = sp_from_expr_with_sum(expr_str, locals_dict, nb_var_expanded)
        else:
            expr_tmp = sp.sympify(expr_str, locals=locals_dict, rational=True)
        
        # Int to float
        expr = expr_tmp
        # expr = sp.N(expr_tmp)
        # Gradient
        gradient = [sp.diff(expr, v) for v in spvars]
        # Hessian (n x n)
        hessian = sp.hessian(expr, spvars)

        list_expr.append([expr]+ list(expr_str1[1:]))
        list_gradient.append([gradient]+ list(expr_str1[1:]))
        list_hessian.append([hessian]+ list(expr_str1[1:]))
        
    # ==========================================================
    # ========================================================== 

    cpp_file = f"{output_file}.hpp"
    path = Path(cpp_file).resolve()
    if not path.exists():
        prepare_output_file(cpp_file,is_gradient_coefficient)

    # Evaluation structure (host and device) written before the class, which uses it.
    # Not extended to gradient expressions: the class then holds its own expressions.
    # Bounds of variables (optional): applied in the evaluation structure
    if bounds and is_gradient_coefficient:
        raise ValueError(f"Error ({class_name}): 'bounds' not extended to gradient expressions.")
    variable_names = [str(v) for v in spvars]
    if has_explicit_variables:
        variable_names += [str(v) for v in explicit_vars]
    if has_auxiliary_variables:
        variable_names += [str(v) for v in auxiliary_vars]
    constant_names = [k for k, v in dict_constants.items() if v != "this->time_"] \
        if has_constants else []
    parsed_bounds = get_bounds(bounds, variable_names, constant_names, class_name)

    eval_generated = False
    if is_gradient_coefficient:
        logging.info(f"{class_name}Eval not generated: not extended to gradient expressions")
    else:
        eval_generated = write_eval_struct(cpp_file, class_name, list_expr, list_gradient,
                                           list_hessian, spvars,
                                           explicit_vars if has_explicit_variables else [],
                                           auxiliary_vars if has_auxiliary_variables else [],
                                           dict_constants if has_constants else {}, coord_vars,
                                           is_conditional_expression, parsed_bounds)

    if eval_generated:
        write_class_using_eval(cpp_file, class_name, list_expr, n, coord_vars,
                               is_conditional_expression)
    else:
        write_class_with_expressions(cpp_file, class_name, list_expr, list_gradient, list_hessian,
                                     spvars, n,
                                     explicit_vars if has_explicit_variables else [],
                                     auxiliary_vars if has_auxiliary_variables else [],
                                     dict_constants if has_constants else {}, coord_vars,
                                     has_explicit_variables, has_auxiliary_variables,
                                     has_constants, has_coordinates, is_gradient_coefficient,
                                     is_conditional_expression)

    logging.info(f"Class {class_name} generated in {cpp_file}")
    # ======================
    # clang-format if found
    # ======================
    clang_format_options = [
        "clang-format",
        "-i",  # Formater le fichier en place
        f"-style={{BasedOnStyle: Google, ColumnLimit: 100, Cpp11BracedListStyle: true, IncludeBlocks: Preserve}}",  
        cpp_file
    ]
    try:
        subprocess.run(clang_format_options, check=True)
        print(f"{cpp_file} is formatted with clang-format")
    except subprocess.CalledProcessError as e:
        print(f"Error when running clang-format: {e}")
    except FileNotFoundError:
        print("clang-format is not found. C++ file not formatted.")
    # ======================


########################################################################################
########################################################################################

def parse_args() -> Namespace:
    """
    Parse and validate command-line arguments for GenerateCoefficient.py.

    This function defines one usage mode: the user provides a JSON file containing 
    coefficient definitions using the `-f/--input-file` option.

    Returns:
        argparse.Namespace: A namespace containing:
            - input_file (Path | None): Path to JSON input file if provided.
            - remove (bool): Whether to remove existing output C++ files.
    """
    def type_file(asstring: str) -> Path:
        path = Path(asstring).resolve()
        if not path.is_file():
            raise ArgumentTypeError(f"file {path} does not exist")
        return path


    parser = ArgumentParser()

    # Group 1: -f is provided
    group1 = parser.add_argument_group('File processing mode')
    group1.add_argument(
        "-f", "--input-file",
        dest="input_file",
        type=type_file,
        help="Path to input file (JSON format, must exist)"
    )

    parser.add_argument("-r", "--remove", dest="remove", help="Remove output cpp files is already exist\n", action="store_true")

    args = parser.parse_args()


    return args


########################################################################################
########################################################################################
########################################################################################

if __name__ == "__main__":
    """
    The script GenerateCoefficient.py performs the following tasks:

    1. Loads coefficient definitions from a JSON file (via `-f`)
    2. Optionally removes previously generated C++ output files (flag `-r`).
    3. Generates C++ classes using the `generate_class_with_functions` function.

    Each coefficient is expected to be a quadruplet:
        (expression, variables, class_name, output_cpp_file)
    """ 
 
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s][%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    args = parse_args()

    
    # Read data
    coefficients = []
    if args.input_file is not None:       
        try:
            with open(args.input_file, 'r') as file:
                data = json.load(file)
            for item in data:
                if "gradient" in item and item["gradient"]:
                    gradient=True
                else:
                    gradient=False
                if "auxiliary_variables" in item:
                    auxiliaries = item["auxiliary_variables"]
                else:
                    auxiliaries = ""
                if "explicit_variables" in item:
                    explicit = item["explicit_variables"]
                else:
                    explicit = ""
                if "constants" in item:
                    constants = item["constants"]
                else:
                    constants = ""
                coordinates = item.get("coordinates", "")
                bounds = item.get("bounds", {})

                expressions = []
                
                cond_var = item.get("range_variable")
                is_conditional_expression = False
                if "expression" in item:
                    expressions.append( (item["expression"], None, None, None, None, None) )
                elif "expressions" in item:
                    is_conditional_expression = True
                    found_default = False
                    if cond_var is None:
                        raise ValueError("Error: 'range_variable' attribute must be defined for a conditional expression.")   
                    if gradient:
                        raise ValueError("Error: conditional expression not extended to gradient expression.")   
                    nb_expression = 0
                    nb_range = 0
                    for expr_item in item["expressions"]:
                        expr = expr_item.get("expression")
                        if expr is not None:
                            nb_expression= nb_expression +1
                        lower = expr_item.get("lower")
                        upper = expr_item.get("upper")
                        lower_strict = expr_item.get("lower_strict")
                        upper_strict = expr_item.get("upper_strict")
                        
                        if expr is None:
                            raise ValueError("Error: expression is expected.")    
                        if not expr_item.get("default"):                  
                            if (lower is not None) and (lower_strict is not None):
                                raise ValueError("Error: either lower or lower_strict must be defined. Please make a choice.")
                            if (lower is None) and (lower_strict is None):
                                raise ValueError("Error: either lower or lower_strict must be defined.")
                            if (upper is not None) and (upper_strict is not None):
                                raise ValueError("Error: either upper or upper_strict must be defined. Please make a choice.")
                            if (upper is None) and (upper_strict is None):
                                raise ValueError("Error: either upper or upper_strict must be defined.")
                            nb_range= nb_range +1
                        else:
                            found_default = True

                        expressions.append( (expr, cond_var, lower, lower_strict, upper, upper_strict) )
                     
                    if not found_default:
                        raise ValueError("Error: a default expression must be defined for each conditional expression .Use the 'default' attribute in the JSON file.")

                    if nb_expression < 2:
                        raise ValueError("Error: at least two expression expressions are expected.")  
                        
                    if nb_range != nb_expression - 1:
                        raise ValueError("Error: the number of ranges must be equal to the number of expressions - 1")

                    check_bounds_consistency(expressions)

                else:
                    raise ValueError("Error: either expression or expressions must be defined.")

                coefficients.append((expressions, item["variables"], explicit, auxiliaries, constants, item["class_name"], item["outputfile"], gradient, is_conditional_expression, coordinates, bounds))
            
        except json.JSONDecodeError:
            print("Error: Failed to decode JSON from the file.")         

    # Remove existing Cpp files
    if args.remove:
        for coef in coefficients:
            output_cpp_file = coef[6]        
            cpp_file = f"{output_cpp_file}.hpp"
            path = Path(cpp_file).resolve()
            if path.is_file():
                os.remove(cpp_file)

    # Generate Cpp files
    for coef in coefficients:
        [expr_str, var_names, explicit_var_names, auxiliary_var_names, constants, class_name, output_cpp_file, is_gradient_coefficient, is_conditional_expression, coordinate_names, bounds] = coef
        generate_class_with_functions(expr_str, var_names, explicit_var_names,auxiliary_var_names, constants, class_name, output_cpp_file, is_gradient_coefficient, is_conditional_expression, coordinate_names, bounds)
        