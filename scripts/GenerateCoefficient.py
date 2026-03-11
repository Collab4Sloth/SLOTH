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
    return dict_constants


#  Case : a * sdot(x(1..j)) = a * (dot(x1,x2)+...+dot(xj,xj))
def sdot_replace_fct(match):
    fact, var, start_str, end_str = match.groups()
    start, end = int(start_str), int(end_str)

    if start > end:
        raise ValueError("Error: invalid bounds in sdot expression. Please check your data.")

    expanded = "+".join(f"dot({var}{i},{var}{i})" for i in range(start, end+1))

    if fact:
        return f"{fact}*({expanded})"
    else:
        return f"{expanded}"
    
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
    f.write(f" * Default expression:\n")
    f.write(f" *     F = {sp.cxxcode(default_expr)}\n")
    f.write(" *\n")

    if is_conditional_expression:
        f.write(" * Conditional definitions:\n")

        for i in range(len(list_expr) - 1):
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

            f.write(f" *   if ({cond_str})\n")
            f.write(f" *       F = {sp.cxxcode(expr)}\n")
            f.write(" *\n")

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
#include <vector>\n
#include "Options/PhysicalPropertiesOptions.hpp"\n  
#include "kernel/Coefficients/FunctionCoefficient.hpp"\n  
#pragma once\n
""")
        
    
def generate_class_with_functions(list_expr_str1, var_names, auxiliary_var_names, constants, class_name, output_file, is_gradient_coefficient, is_conditional_expression):
    
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

    with open(cpp_file, "a") as f:

        f.write("/**\n")
        f.write(" *\n")
        f.write(f" * @brief C++ function of the expression\n")
        f.write(" *\n")

        print_header(f, list_expr, is_conditional_expression)
 
        # Class
        f.write(f"class {class_name} : public FunctionCoefficient {{\n")
        f.write(" private:\n")
        f.write("  double prefactor_;\n")
        f.write(" protected:\n")
        f.write("  std::function<double(const std::span<const double>&,const std::span<const double>&, const unsigned int dimension)> F() final;\n")
        f.write("  std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const unsigned int dimension)> GradientF() final;\n")
        f.write("  std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const unsigned int dimension)> HessianF() final;\n\n")
        f.write(" public:\n")
        f.write(f"  {class_name}() {{this->prefactor_ = 1.0; }}\n")
        f.write(f" explicit {class_name}(const double prefactor) {{this->prefactor_ = prefactor; }}\n")
        f.write(f"  virtual ~{class_name}() = default;\n")
        f.write("};\n\n")


        f.write("/**\n")
        f.write(" *\n")
        f.write(f" * @brief C++ function of the expression\n")
        f.write(" *\n")

        # Expression par défaut
        default_expr = list_expr[-1][0]
        f.write(f" *   F = {sp.cxxcode(default_expr)}\n")
        f.write(f" *\n")
        for i in range(0, len(list_expr)-1):
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

            f.write(f" *   if ({cond_str})\n")
            f.write(f" *       F = {sp.cxxcode(expr)}\n")
            f.write(f" * \n")
        f.write(f" * @return std::function<double(const std::span<const double>&,const std::span<const double>&)>\n")
        f.write(f" */\n ")  
        
        # Fonction F()
        f.write(f"std::function<double(const std::span<const double>&,const std::span<const double>&, const unsigned int dimension)> {class_name}::F() {{\n")

        if(not is_gradient_coefficient):        
            if(has_auxiliary_variables):
                f.write("  auto func = [&](const std::span<const double>& input_vector, const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {\n")
            else:
                f.write("  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const unsigned int dimension) {\n")

            for i, v in enumerate(spvars):
                f.write(f"    double {v} = input_vector[{i}];\n")
            if(has_auxiliary_variables):
                for i, v in enumerate(auxiliary_vars):
                    f.write(f"    double {v} = auxiliary_vector[{i}];\n")
            
            if(has_constants):
                for const in dict_constants:
                    f.write(f"    double {const} = {dict_constants[const]};\n")
   
            #====================================================
            #====================================================
            # single expression or the default one in case "conditional expression"
            expr_0 = list_expr[-1][0]
            f.write(f"    double F = {sp.cxxcode(expr_0)};\n")
            if is_conditional_expression:
                for i in range(0,len(list_expr)-1):
                    expr = list_expr[i][0]
                    # ===============
                    # <= && >=
                    if list_expr[i][2] is not None and list_expr[i][4] is not None:
                        f.write(f"   if(({list_expr[i][2]}<={list_expr[i][1]}) && ({list_expr[i][1]}<={list_expr[i][4]}))  \n")
                        f.write("{ \n")
                        f.write(f"    F = {sp.cxxcode(expr)};\n")
                        f.write("} \n")
                    # ===============
                    # < && >
                    if list_expr[i][3] is not None and list_expr[i][5] is not None:
                        f.write(f"   if(({list_expr[i][3]}<{list_expr[i][1]}) && ({list_expr[i][1]}<{list_expr[i][5]}))  \n")
                        f.write("{ \n")
                        f.write(f"    F = {sp.cxxcode(expr)};\n")
                        f.write("} \n")
                    # ===============
                    # <= && >
                    if list_expr[i][2] is not None and list_expr[i][5] is not None:
                        f.write(f"   if(({list_expr[i][2]}<={list_expr[i][1]}) && ({list_expr[i][1]}<{list_expr[i][5]}))  \n")
                        f.write("{ \n")
                        f.write(f"    F = {sp.cxxcode(expr)};\n")
                        f.write("} \n")
                    # ===============
                    # < && >=
                    if list_expr[i][3] is not None and list_expr[i][4] is not None:
                        f.write(f"   if(({list_expr[i][3]}<{list_expr[i][1]}) && ({list_expr[i][1]}<={list_expr[i][4]}))  \n")
                        f.write("{ \n")
                        f.write(f"    F = {sp.cxxcode(expr)};\n")
                        f.write("} \n")
                
        else:
            if(has_auxiliary_variables):
                f.write("  auto func = [&](const std::span<const double>& input_vector, const std::span<const double>& auxiliary_vector, const unsigned int dimension) {\n")
            else:
                f.write("  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&, const unsigned int dimension) {\n")

            for i, v in enumerate(spvars):
                f.write(f"    std::vector<double> {v};\n")
                f.write(f"    for(unsigned int i=0;i<dimension;i++) {v}.push_back(input_vector[{i}*dimension+i]);\n")

            if(has_auxiliary_variables):
                for i, v in enumerate(auxiliary_vars):
                    f.write(f"    std::vector<double> {v};\n")
                    f.write(f"    for(unsigned int i=0;i<dimension;i++) {v}.push_back(auxiliary_vector[{i}*dimension+i]);\n")

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
        f.write(f"std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const unsigned int dimension)> {class_name}::GradientF() {{\n")


        if(not is_gradient_coefficient):   
            if(has_auxiliary_variables):
                f.write("  auto func = [&](const std::span<const double>& input_vector, const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {\n")
            else:
                f.write("  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const unsigned int dimension) {\n")

            for i, v in enumerate(spvars):
                f.write(f"    double {v} = input_vector[{i}];\n")
            if(has_auxiliary_variables):
                for i, v in enumerate(auxiliary_vars):
                    f.write(f"    double {v} = auxiliary_vector[{i}];\n")
                    
            if(has_constants):
                for const in dict_constants:
                    f.write(f"    double {const} = {dict_constants[const]};\n")
                    
            f.write(f"    std::vector<double> gradient({n});\n")

            #====================================================
            #====================================================
            # single gradient or the default one in case "conditional expression"
            gradient_0 = list_gradient[-1][0]
            print_gradient(f, gradient_0, n)

            if is_conditional_expression:
                for i in range(0,len(list_gradient)-1):
                    gradient = list_gradient[i][0]
                    # ===============
                    # <= && >=
                    if list_gradient[i][2] is not None and list_gradient[i][4] is not None:
                        f.write(f"   if(({list_gradient[i][2]}<={list_gradient[i][1]}) && ({list_gradient[i][1]}<={list_gradient[i][4]}))  \n")
                        f.write("{ \n")
                        print_gradient(f, gradient, n)
                        f.write("} \n")
                    # ===============
                    # < && >
                    if list_gradient[i][3] is not None and list_gradient[i][5] is not None:
                        f.write(f"   if(({list_gradient[i][3]}<{list_gradient[i][1]}) && ({list_gradient[i][1]}<{list_gradient[i][5]}))  \n")
                        f.write("{ \n")
                        print_gradient(f, gradient, n)
                        f.write("} \n")
                    # ===============
                    # <= && >
                    if list_gradient[i][2] is not None and list_gradient[i][5] is not None:
                        f.write(f"   if(({list_gradient[i][2]}<={list_gradient[i][1]}) && ({list_gradient[i][1]}<{list_gradient[i][5]}))  \n")
                        f.write("{ \n")
                        print_gradient(f, gradient, n)
                        f.write("} \n")
                    # ===============
                    # < && >=
                    if list_gradient[i][3] is not None and list_gradient[i][4] is not None:
                        f.write(f"   if(({list_gradient[i][3]}<{list_gradient[i][1]}) && ({list_gradient[i][1]}<={list_gradient[i][4]}))  \n")
                        f.write("{ \n")
                        print_gradient(f, gradient, n)
                        f.write("} \n")
        else: 
            f.write("  auto func = [&]([[maybe_unused]] const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const unsigned int dimension) {\n")
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
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const unsigned int dimension)> 
 */
""")
        # HessianF()
        f.write(f"std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const unsigned int dimension)> {class_name}::HessianF() {{\n")
        if(not is_gradient_coefficient): 
            if(has_auxiliary_variables):
                f.write("  auto func = [&](const std::span<const double>& input_vector, const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {\n")
            else:
                f.write("  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const unsigned int dimension) {\n")
            for i, v in enumerate(spvars):
                f.write(f"    double {v} = input_vector[{i}];\n")
            if(has_auxiliary_variables):
                for i, v in enumerate(auxiliary_vars):
                    f.write(f"    double {v} = auxiliary_vector[{i}];\n")
                
            if(has_constants):
                for const in dict_constants:
                    f.write(f"    double {const} = {dict_constants[const]};\n")
                    
            f.write(f"    std::vector<double> hessian({n*n});\n")
       
            #====================================================
            #====================================================
            # single hessian or the default one in case "conditional expression"
            hessian_0 = list_hessian[-1][0]
            print_hessian(f, hessian_0, n)
                        
            if is_conditional_expression:
                for i in range(0,len(list_hessian)-1):
                    hessian = list_hessian[i][0]
                    # ===============
                    # <= && >=
                    if list_hessian[i][2] is not None and list_hessian[i][4] is not None:
                        f.write(f"   if(({list_hessian[i][2]}<={list_hessian[i][1]}) && ({list_hessian[i][1]}<={list_hessian[i][4]}))  \n")
                        f.write("{ \n")
                        print_hessian(f, hessian, n)
                        f.write("} \n")
                    # ===============
                    # < && >
                    if list_hessian[i][3] is not None and list_hessian[i][5] is not None:
                        f.write(f"   if(({list_hessian[i][3]}<{list_hessian[i][1]}) && ({list_hessian[i][1]}<{list_hessian[i][5]}))  \n")
                        f.write("{ \n")
                        print_hessian(f, hessian, n)
                        f.write("} \n")
                    # ===============
                    # <= && >
                    if list_hessian[i][2] is not None and list_hessian[i][5] is not None:
                        f.write(f"   if(({list_hessian[i][2]}<={list_hessian[i][1]}) && ({list_hessian[i][1]}<{list_hessian[i][5]}))  \n")
                        f.write("{ \n")
                        print_hessian(f, hessian, n)
                        f.write("} \n")
                    # ===============
                    # < && >=
                    if list_hessian[i][3] is not None and list_hessian[i][4] is not None:
                        f.write(f"   if(({list_hessian[i][3]}<{list_hessian[i][1]}) && ({list_hessian[i][1]}<={list_hessian[i][4]}))  \n")
                        f.write("{ \n")
                        print_hessian(f, hessian, n)
                        f.write("} \n")

        else:
            f.write("  auto func = [&]([[maybe_unused]] const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&, [[maybe_unused]] const unsigned int dimension) {\n")
            f.write(f"    std::vector<double> hessian({n*n},0.0);\n")

        f.write("    return hessian;\n")
        f.write("  };\n")
        f.write("  return func;\n")
        f.write("}\n")


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
                if "constants" in item:
                    constants = item["constants"]
                else:
                    constants = ""

                expressions = []
                is_conditional_expression = False
                if "expression" in item:
                    expressions.append( (item["expression"], None, None, None, None, None) )
                elif "expressions" in item:
                    is_conditional_expression = True
                    if gradient:
                        raise ValueError("Error: conditional expression not extended to gradient expression.")   
                    nb_expression = 0
                    nb_condition = 0
                    for expr_item in item["expressions"]:
                        expr = expr_item.get("expression")
                        if expr is not None:
                            nb_expression= nb_expression +1
                        cond_var = expr_item.get("range_variable")
                        if cond_var is not None:
                            nb_condition= nb_condition +1
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

                        expressions.append( (expr, cond_var, lower, lower_strict, upper, upper_strict) )
                        
                    if nb_expression < 2:
                        raise ValueError("Error: at least two expression expressions are expected.")  
                        
                    if nb_condition < 1:
                        raise ValueError("Error: at least one condition is expected.")       

                else:
                    raise ValueError("Error: either expression or expressions must be defined.")

                coefficients.append((expressions, item["variables"], auxiliaries, constants, item["class_name"], item["outputfile"], gradient, is_conditional_expression))
            
        except json.JSONDecodeError:
            print("Error: Failed to decode JSON from the file.")         

    # Remove existing Cpp files
    if args.remove:
        for coef in coefficients:
            output_cpp_file = coef[-3]        
            cpp_file = f"{output_cpp_file}.hpp"
            path = Path(cpp_file).resolve()
            if path.is_file():
                os.remove(cpp_file)

    # Generate Cpp files
    for coef in coefficients:
        [expr_str, var_names, auxiliary_var_names, constants, class_name, output_cpp_file, is_gradient_coefficient, is_conditional_expression] = coef
        generate_class_with_functions(expr_str, var_names, auxiliary_var_names, constants, class_name, output_cpp_file, is_gradient_coefficient, is_conditional_expression)
    
