#!/usr/bin/env python3
"""
List the SLOTH kernel functions never called, from a gcovr JSON report.

Usage:
    uncovered_functions.py coverage.json                 # summary per module
    uncovered_functions.py coverage.json -m Parameters   # details for one module
    uncovered_functions.py coverage.json --all           # details for all modules

The JSON report is produced by gcovr --json (see the unit_coverage target in
tests/Unit/CMakeLists.txt).

A template function has one entry per instantiation in gcov (one per DIM, per
type...). Entries are grouped by (file, line): a function is "covered" if at
least one of its instantiations has been called.
"""
import argparse
import json
import sys
from collections import defaultdict

KERNEL_DIRS = ("kernel/include/", "kernel/src/")


def module_of(path):
    """kernel/include/Parameters/Parameters.tpp -> ('Parameters', 'Parameters/Parameters.tpp')"""
    for prefix in KERNEL_DIRS:
        pos = path.find(prefix)
        if pos >= 0:
            rel = path[pos + len(prefix):]
            return rel.split("/")[0], rel
    return None, path


def function_name(entry):
    # gcovr >= 7: "demangled_name"; older versions: "name"
    return entry.get("demangled_name") or entry.get("name") or "?"


def function_count(entry):
    # gcovr >= 7: "execution_count"; older versions: "count"
    return entry.get("execution_count", entry.get("count", 0)) or 0


def collect(report):
    """{module: {(rel_path, line): [name, total_count]}}"""
    modules = defaultdict(dict)
    for f in report.get("files", []):
        module, rel = module_of(f.get("file", ""))
        if module is None:
            continue
        for fn in f.get("functions", []):
            key = (rel, fn.get("lineno", 0))
            name = function_name(fn)
            count = function_count(fn)
            if key in modules[module]:
                modules[module][key][1] += count
            else:
                modules[module][key] = [name, count]
    return modules


def print_summary(modules):
    rows = []
    for module, funcs in modules.items():
        total = len(funcs)
        uncovered = sum(1 for _, count in funcs.values() if count == 0)
        rows.append((module, total, uncovered))
    # Modules with the most uncovered functions first
    rows.sort(key=lambda r: (-r[2], r[0]))

    print(f"{'Module':<24}{'functions':>10}{'uncovered':>11}{'covered %':>11}")
    print("-" * 56)
    all_total = all_uncovered = 0
    for module, total, uncovered in rows:
        pct = 100.0 * (total - uncovered) / total if total else 100.0
        print(f"{module:<24}{total:>10}{uncovered:>11}{pct:>10.1f}%")
        all_total += total
        all_uncovered += uncovered
    print("-" * 56)
    pct = 100.0 * (all_total - all_uncovered) / all_total if all_total else 100.0
    print(f"{'TOTAL':<24}{all_total:>10}{all_uncovered:>11}{pct:>10.1f}%")


def print_details(modules, selected):
    for module in sorted(selected):
        funcs = modules.get(module)
        if funcs is None:
            print(f"Unknown module '{module}'. Available: {', '.join(sorted(modules))}")
            continue
        uncovered = sorted((key, name) for key, (name, count) in funcs.items() if count == 0)
        print(f"\n=== {module}: {len(uncovered)} uncovered function(s) out of {len(funcs)}")
        for (path, line), name in uncovered:
            print(f"  {path}:{line}\n      {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("json_report", help="gcovr JSON report (gcovr --json)")
    parser.add_argument("-m", "--module", action="append", default=[],
                        help="details for this module (repeatable)")
    parser.add_argument("--all", action="store_true", help="details for all modules")
    args = parser.parse_args()

    with open(args.json_report, encoding="utf-8") as stream:
        report = json.load(stream)
    modules = collect(report)
    if not modules:
        sys.exit("No kernel function found in the report (check the gcovr --filter options)")

    print_summary(modules)
    selected = set(modules) if args.all else set(args.module)
    if selected:
        print_details(modules, selected)


if __name__ == "__main__":
    main()
