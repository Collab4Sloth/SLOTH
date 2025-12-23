#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

    


def prepare_output_file(output_file, is_gradient_coefficient):

    filename = Path(output_file).resolve()

    with open(output_file, "w") as f:
        # Header 
        f.write("""/**
 *
 * Copyright CEA (C) 2025
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
#include <cmath>
#include <functional>
#include <numeric>
#include <vector>\n
#include "kernel/Coefficients/FunctionCoefficient.hpp"\n  
#pragma once\n
""")

def dot_to_inner_product(code: str) -> str:
    """
    Replace Sum(a[i]*b[i], (i, 0, n)) by std::inner_product(a.begin(), a.end(), b.begin(), 0)
    """
    pattern = re.compile(
        r"dot\(\s*"        # dot(
        r"(\w+)\s*,\s*"    # première variable
        r"(\w+)\s*\)"      # deuxième variable
    )

    def repl(match):
        var1 = match.group(1)
        var2 = match.group(2)
        return f"std::inner_product({var1}.begin(), {var1}.end(), {var2}.begin(), 0.0)"

    return pattern.sub(repl, code)

def generate_class_with_functions(expr_str, var_names, auxiliary_var_names, class_name, output_file, is_gradient_coefficient):
    var_names+=","
    vars = sp.symbols(var_names)
    n = len(vars)
    has_auxiliary_variables = False
    if(len(auxiliary_var_names)>0):
        has_auxiliary_variables = True
        auxiliary_var_names+=","
        auxiliary_vars = sp.symbols(auxiliary_var_names)
        print(auxiliary_vars)

    locals_dict={}
    if(is_gradient_coefficient):
        dot = sp.Function('dot')
        locals_dict["dot"]=dot

    expr_tmp = sp.sympify(expr_str, locals=locals_dict, rational=True)


    # Int to float
    expr = expr_tmp
    # expr = sp.N(expr_tmp)
    # Gradient
    gradient = [ sp.diff(expr, v) for v in vars]
    # Hessian (n x n)
    hessian = sp.hessian(expr, vars)

    cpp_file = f"{output_file}.hpp"
    path = Path(cpp_file).resolve()
    if not path.exists():
        prepare_output_file(cpp_file,is_gradient_coefficient)

    with open(cpp_file, "a") as f:

        f.write(f"""
/**
 *
 * @brief Coefficient based on expression: {expr_str}
 *
 */
""")
        # Class
        f.write(f"class {class_name} : public FunctionCoefficient {{\n")
        f.write(" private:\n")
        f.write("  double prefactor_;\n")
        f.write(" protected:\n")
        f.write("  std::function<double(const std::vector<double>&,const std::vector<double>&, const int dimension)> F() final;\n")
        f.write("  std::function<std::vector<double>(const std::vector<double>&,const std::vector<double>&, const int dimension)> GradientF() final;\n")
        f.write("  std::function<std::vector<double>(const std::vector<double>&,const std::vector<double>&, const int dimension)> HessianF() final;\n\n")
        f.write(" public:\n")
        f.write(f"  {class_name}() {{this->prefactor_ = 1.0; }};\n")
        f.write(f"  {class_name}(const double prefactor) {{this->prefactor_ = prefactor; }};\n")
        f.write(f"  ~{class_name}(){{}};\n")
        f.write("};\n\n")


        f.write(f"""/**
 *
 * @brief C++ function of the expression: {expr_str}
 * 
 * @return std::function<double(const std::vector<double>&,const std::vector<double>&)> 
 */
""")
        # Fonction F()
        f.write(f"std::function<double(const std::vector<double>&,const std::vector<double>&, const int dimension)> {class_name}::F() {{\n")

        if(not is_gradient_coefficient):        
            if(has_auxiliary_variables):
                f.write("  auto func = [&](const std::vector<double>& input_vector, const std::vector<double>& auxiliary_vector, [[maybe_unused]] const int dimension) {\n")
            else:
                f.write("  auto func = [&](const std::vector<double>& input_vector, [[maybe_unused]] const std::vector<double>&, [[maybe_unused]] const int dimension) {\n")

            for i, v in enumerate(vars):
                f.write(f"    double {v} = input_vector[{i}];\n")
            if(has_auxiliary_variables):
                for i, v in enumerate(auxiliary_vars):
                    f.write(f"    double {v} = auxiliary_vector[{i}];\n")
            f.write(f"    double F = {sp.cxxcode(expr)};\n")
        else:
            if(has_auxiliary_variables):
                f.write("  auto func = [&](const std::vector<double>& input_vector, const std::vector<double>& auxiliary_vector, const int dimension) {\n")
            else:
                f.write("  auto func = [&](const std::vector<double>& input_vector, [[maybe_unused]] const std::vector<double>&, const int dimension) {\n")

            for i, v in enumerate(vars):
                f.write(f"    std::vector<double> {v};\n")
                f.write(f"    for(unsigned int i=0;i<dimension;i++) {v}.push_back(input_vector[{i}*dimension+i]);\n")

            if(has_auxiliary_variables):
                for i, v in enumerate(auxiliary_vars):
                    f.write(f"    std::vector<double> {v};\n")
                    f.write(f"    for(unsigned int i=0;i<dimension;i++) {v}.push_back(auxiliary_vector[{i}*dimension+i]);\n")


            f.write(f"    double F = {dot_to_inner_product(str(sp.N(expr)))};\n")

   
        f.write("    return this->prefactor_ * F;\n")
        f.write("  };\n")
        f.write("  return func;\n")
        f.write("}\n\n")

        f.write(f"""/**
 *
 * @brief Gradient
 * 
 * @return std::function<std::vector<double>(const std::vector<double>&,const std::vector<double>&, const int dimension)> 
 */
""")
        # Gradient
        f.write(f"std::function<std::vector<double>(const std::vector<double>&,const std::vector<double>&, const int dimension)> {class_name}::GradientF() {{\n")


        if(not is_gradient_coefficient):   
            if(has_auxiliary_variables):
                f.write("  auto func = [&](const std::vector<double>& input_vector, const std::vector<double>& auxiliary_vector, [[maybe_unused]] const int dimension) {\n")
            else:
                f.write("  auto func = [&](const std::vector<double>& input_vector, [[maybe_unused]] const std::vector<double>&, [[maybe_unused]] const int dimension) {\n")

            for i, v in enumerate(vars):
                f.write(f"    double {v} = input_vector[{i}];\n")
            if(has_auxiliary_variables):
                for i, v in enumerate(auxiliary_vars):
                    f.write(f"    double {v} = auxiliary_vector[{i}];\n")
            f.write(f"    std::vector<double> gradient({n});\n")
            for i in range(n):
                f.write(f"    gradient[{i}] = this->prefactor_ * ({sp.cxxcode(gradient[i])});\n")
        else: 
            f.write("  auto func = [&]([[maybe_unused]] const std::vector<double>& input_vector, [[maybe_unused]] const std::vector<double>&, [[maybe_unused]] const int dimension) {\n")
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
 * @return std::function<std::vector<double>(const std::vector<double>&,const std::vector<double>&, const int dimension)> 
 */
""")
        # HessianF()
        f.write(f"std::function<std::vector<double>(const std::vector<double>&,const std::vector<double>&, const int dimension)> {class_name}::HessianF() {{\n")
        if(not is_gradient_coefficient): 
            if(has_auxiliary_variables):
                f.write("  auto func = [&](const std::vector<double>& input_vector, const std::vector<double>& auxiliary_vector, [[maybe_unused]] const int dimension) {\n")
            else:
                f.write("  auto func = [&](const std::vector<double>& input_vector, [[maybe_unused]] const std::vector<double>&, [[maybe_unused]] const int dimension) {\n")
            for i, v in enumerate(vars):
                f.write(f"    double {v} = input_vector[{i}];\n")
            if(has_auxiliary_variables):
                for i, v in enumerate(auxiliary_vars):
                    f.write(f"    double {v} = auxiliary_vector[{i}];\n")
                
            f.write(f"    std::vector<double> hessian({n*n});\n")
            for i in range(n):
                for j in range(n):
                    f.write(f"    hessian[{i*n + j}] = this->prefactor_ * ({sp.cxxcode(hessian[i,j])});\n")
        else:
            f.write("  auto func = [&]([[maybe_unused]] const std::vector<double>& input_vector, [[maybe_unused]] const std::vector<double>&, [[maybe_unused]] const int dimension) {\n")
            f.write(f"    std::vector<double> hessian({n*n},0.0);\n")


        f.write("    return hessian;\n")
        f.write("  };\n")
        f.write("  return func;\n")
        f.write("}\n")

    print(f"Class {class_name} generated in {cpp_file}")

########################################################################################
########################################################################################

def parse_args() -> Namespace:
    """
    Parse and validate command-line arguments for GenerateCoefficient.py.

    This function defines two mutually exclusive usage modes:
    
    1. **File processing mode**:  
       The user provides a JSON file containing coefficient definitions using
       the `-f/--input-file` option.

    2. **Direct processing mode**:  
       The user provides one or several coefficient definitions directly from
       the command line using `-c/--coefficients`.  
       Each coefficient must follow the format: ``expr,var,class_name,outputfile``.

    Several validations are performed:
      - ``-f`` must point to an existing file.
      - ``-c`` must contain valid comma-separated quadruplets.
      - ``-c`` is required if ``-f`` is omitted.

    Returns:
        argparse.Namespace: A namespace containing:
            - input_file (Path | None): Path to JSON input file if provided.
            - coefficients (list[tuple[str, str, str, str]] | None):
              Coefficients provided directly via command line.
            - remove (bool): Whether to remove existing output C++ files.
    """
    def type_file(asstring: str) -> Path:
        path = Path(asstring).resolve()
        if not path.is_file():
            raise ArgumentTypeError(f"file {path} does not exist")
        return path

    def type_coefficient(val: str) -> Tuple[str, str, str, str]:
        """Convert string in format 'a,b,c,d' to tuple (a,b,c,d)"""
        parts = [p.strip() for p in val.split(',')]
        if len(parts) != 4:
            raise ArgumentTypeError(f"Expected format 'a,b,c,d', got '{val}'")
        return tuple(parts)

    parser = ArgumentParser()

    # Group 1: -f is provided (optional -c)
    group1 = parser.add_argument_group('File processing mode')
    group1.add_argument(
        "-f", "--input-file",
        dest="input_file",
        type=type_file,
        help="Path to input file (JSON format, must exist)"
    )

    # Group 2: -f is not provided (required -c)
    group2 = parser.add_argument_group('Direct processing mode')
    group2.add_argument(
        "-c", "--coefficients",
        dest="coefficients",
        nargs='+',
        type=type_coefficient,
        help="Space-separated list of strings (e.g., 'expr1,var1,class1,cppfile1' 'expr2,var2,class2,cppfile2')"
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

    1. Loads coefficient definitions either from:
    - a JSON file (via `-f`), or
    - command-line tuples (via `-c`).
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
                
                coefficients.append((item["expression"], item["variables"], auxiliaries, item["class_name"], item["outputfile"], gradient))
            


        except json.JSONDecodeError:
            print("Error: Failed to decode JSON from the file.")           
    else:
        # Direct mode: -f is not provided, -c is required
        if args.coefficients is None:
            parser.error("When -f is not provided, -c is required")

        coefficients = args.coefficients

    # Remove existing Cpp files
    if args.remove:
        for coef in coefficients:
            output_cpp_file = coef[-2]        
            cpp_file = f"{output_cpp_file}.hpp"
            path = Path(cpp_file).resolve()
            if path.is_file():
                os.remove(cpp_file)

    # Generate Cpp files
    for coef in coefficients:
        [expr_str, var_names, auxiliary_var_names, class_name, output_cpp_file, is_gradient_coefficient] = coef
        generate_class_with_functions(expr_str, var_names, auxiliary_var_names, class_name, output_cpp_file, is_gradient_coefficient)
    
