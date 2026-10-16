#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare columns of two CSV files (results vs reference) with an error criterion.

Comparison works with any Python 3 (system, Homebrew, Spack, venv, conda).
"""

import csv
import math
import sys
from argparse import ArgumentParser, RawTextHelpFormatter


def calculate_error(a, b, error_type):
    if error_type == 'relativeEpsilon':
        epsilon = 1.e-20
        return _divide(abs(b - a), abs(max(a, epsilon)))
    elif error_type == 'relative':
        return _divide(abs(b - a), abs(a))
    elif error_type == 'absolute':
        return abs(b - a)
    else:
        raise ValueError(
            "Error type must be either 'relative', 'relativeEpsilon', or 'absolute'")


def _divide(num, den):
    # x/0 -> inf (difference detected), 0/0 -> nan (no difference detected)
    if den == 0:
        return math.nan if (num == 0 or math.isnan(num)) else math.inf
    return num / den


def _to_float(value):
    # Empty cell -> nan (missing value, never reported as a difference)
    value = value.strip()
    if value == '':
        return math.nan
    return float(value)


def read_csv(filename):
    """Return (headers, rows) with rows as lists of floats; blank lines are skipped."""
    with open(filename, newline='') as f:
        lines = [row for row in csv.reader(f) if any(cell.strip() for cell in row)]
    if not lines:
        raise ValueError(f"Error : empty file '{filename}'")
    headers = lines[0]
    rows = [[_to_float(cell) for cell in row] for row in lines[1:]]
    return headers, rows


def compare_csv_files(file1, file2, columns, error_type, criterion, skip_lines):
    headers1, rows1 = read_csv(file1)
    _, rows2 = read_csv(file2)

    # Full comparison if a negative value is provided
    if columns[0] < 0:
        columns = list(range(len(headers1)))

    errors = {i + skip_lines: {} for i in range(len(rows1))}
    nrows = min(len(rows1), len(rows2))

    # The first listed column (the iteration/time column) is not compared
    for id, col in enumerate(columns):
        col_name = headers1[col]
        for i in range(nrows):
            if id > 0:
                error = calculate_error(rows1[i][col], rows2[i][col], error_type)
                if error > criterion:
                    errors[i + skip_lines][col_name] = error

    errors = {line: cols for line, cols in errors.items() if cols}

    if errors:
        generate_markdown_report(errors, file1, file2, error_type, criterion)
        raise ValueError(
            "Error : differences are detected, see comparison_report file for details")


def generate_markdown_report(errors, file1, file2, error_type, criterion):
    report_file = 'comparison_report.md'
    with open(report_file, 'w') as f:
        f.write("# Comparison Report\n\n")
        f.write(f"Comparing files: `{file1}` and `{file2}`\n\n")
        f.write(f"Error type: `{error_type}`\n")
        f.write(f"Error criterion: `{criterion}`\n\n")
        f.write("## Detected Differences\n\n")

        columns = sorted({col for cols in errors.values() for col in cols})
        header = "| Iter[-] | " + " | ".join(columns) + " |\n"
        separator = "|-------------|" + "-------|" * len(columns) + "\n"

        f.write(header)
        f.write(separator)

        for line_number, cols in errors.items():
            row = [str(line_number)]
            for col in columns:
                row.append(str(cols.get(col, "")))
            f.write("| " + " | ".join(row) + " |\n")

    print(f"Markdown report generated: {report_file}")


if __name__ == '__main__':
    parser = ArgumentParser(description='Compare columns of two CSV files with specified error criterion',
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Increase output verbosity")
    parser.add_argument("-f", "--files", nargs=2,
                        required=True, help="List of files to compare : reference, prediction")
    parser.add_argument("-c", "--columns", nargs='*', type=int, required=True,
                        help="List of columns to compare (0-indexed, negative value for all)")
    parser.add_argument("-e", "--error", required=True, choices=['relative', 'relativeEpsilon', 'absolute'],
                        help="Type of error: 'relative', 'relativeEpsilon', or 'absolute'")
    parser.add_argument("-k", "--criterion", type=float,
                        required=True, help="Criterion to compare data")
    parser.add_argument("-w", "--switch", type=int,
                        default=0, help="Number of lines to ignore")

    args = parser.parse_args()

    if args.verbose:
        print(f"Comparing files: {args.files[0]} and {args.files[1]}")
        print(f"Columns to compare: {args.columns}")
        print(f"Error type: {args.error}")
        print(f"Error criterion: {args.criterion}")
        print(f"Lines to skip: {args.switch}")

    try:
        compare_csv_files(
            args.files[0], args.files[1], args.columns, args.error, args.criterion, args.switch)
    except ValueError as e:
        print(e)
        sys.exit(1)

    if args.verbose:
        print("Comparison successful")
