#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter


def calculate_error(a, b, error_type):
    if error_type == 'relativeEpsilon':
        epsilon = 1.e-20
        return abs((a - b) / max(a, b, epsilon))
    elif error_type == 'relative':
        return abs((a - b) / max(a, b))
    elif error_type == 'absolute':
        return abs(a - b)
    else:
        raise ValueError(
            "Error type must be either 'relative', 'relativeEpsilon', or 'absolute'")


def compare_csv_files(file1, file2, columns, error_type, criterion, skip_lines):
    # Load CSV files into DataFrames
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    # Headers
    headers1 = df1.columns.tolist()
    headers2 = df2.columns.tolist()

    # Compare headers
    header_errors = []
    if headers1 != headers2:
        for i, (h1, h2) in enumerate(zip(headers1, headers2)):
            if h1 != h2:
                header_errors.append((i, h1, h2))

    # Full comparison if a negative value is provided
    if columns[0] < 0:
        columns = list(range(len(df1.columns)))

    errors = {i + skip_lines: {} for i in range(len(df1))}
    generate_md = False
    for id, col in enumerate(columns):
        col_name = headers1[col]
        for i in range(len(df1)):
            if id > 0:
                error = calculate_error(
                    df1.iloc[i, col], df2.iloc[i, col], error_type)
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
        f.write(f"# Comparison Report\n\n")
        f.write(f"Comparing files: `{file1}` and `{file2}`\n\n")
        f.write(f"Error type: `{error_type}`\n")
        f.write(f"Error criterion: `{criterion}`\n\n")
        f.write(f"## Detected Differences\n\n")

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
                        required=True, help="List of files to compare")
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
