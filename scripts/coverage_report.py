#!/usr/bin/env python3
"""
Coverage report of the SLOTH kernel (target `make cc`).

Usage:
    coverage_report.py --build BUILD_DIR --source SOURCE_DIR --output OUTPUT_DIR
                       [--gcov-executable "xcrun llvm-cov gcov"]

Reads the coverage data of the tests that have been run (ctest), and writes in
OUTPUT_DIR:
    coverage.html (+ one page per file)  detailed report
    coverage.xml                         Cobertura format
    coverage.json                        input of uncovered_functions.py
then prints the number of uncovered kernel functions per module. On GitHub
Actions, the same table is added to the job summary ($GITHUB_STEP_SUMMARY).

Only the object directories (CMakeFiles/<target>.dir) containing coverage
counters (*.gcda) are read: the targets that have been built but not run, or
built from older versions of the sources (renamed tests, tests not rebuilt
after a change of a kernel header), are ignored. In these directories, the
objects without counters (source files removed from their target) are ignored
too. The code of kernel/include
compiled in the test executables (templates, inline functions) is taken into
account.
"""
import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

import uncovered_functions


def object_dirs_with_counters(build_dir):
    """CMakeFiles/<target>.dir directories of the build tree containing *.gcda files"""
    dirs = set()
    for root, _, files in os.walk(build_dir):
        if not any(name.endswith(".gcda") for name in files):
            continue
        path = Path(root)
        # Object directory of the target: the CMakeFiles/<target>.dir ancestor
        for parent in [path, *path.parents]:
            if parent.suffix == ".dir" and parent.parent.name == "CMakeFiles":
                dirs.add(parent)
                break
    return sorted(dirs)


def orphan_notes(dirs):
    """*.gcno files without *.gcda in the selected object directories.

    When an executable runs, a *.gcda file is written for every object it
    contains: a *.gcno file alone is an object that is no longer part of the
    executable (e.g. a source file removed from the target).
    """
    orphans = []
    for directory in dirs:
        for gcno in directory.rglob("*.gcno"):
            if not gcno.with_suffix(".gcda").exists():
                orphans.append(gcno)
    return orphans


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--build", required=True, help="build directory")
    parser.add_argument("--source", required=True, help="SLOTH source directory")
    parser.add_argument("--output", required=True, help="output directory of the reports")
    parser.add_argument("--gcov-executable", default=None,
                        help="gcov command matching the compiler (e.g. 'gcov-14')")
    args = parser.parse_args()

    source = Path(args.source).resolve()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)

    dirs = object_dirs_with_counters(Path(args.build).resolve())
    if not dirs:
        sys.exit("No coverage counters (*.gcda) found: run the tests first (ctest)")
    print(f"Coverage counters found in {len(dirs)} target(s)")
    orphans = orphan_notes(dirs)
    for gcno in orphans:
        print(f"Ignored (object no longer part of its target): {gcno}")

    command = ["gcovr", "-r", str(source), *map(str, dirs),
               "--filter", str(source / "kernel" / "src"),
               "--filter", str(source / "kernel" / "include"),
               "--exclude", str(source / "kernel" / "include" / "MAToolsProfiling") + "/.*",
               "--exclude-unreachable-branches",
               "--json", str(output / "coverage.json"),
               "--xml", str(output / "coverage.xml"),
               "--html-details", str(output / "coverage.html"),
               "--print-summary",
               # Counters of a source file removed since they were written
               "--gcov-ignore-errors", "source_not_found"]
    for gcno in orphans:
        command += ["--gcov-exclude", re.escape(str(gcno))]
    if args.gcov_executable:
        command += ["--gcov-executable", args.gcov_executable]
    print(" ".join(shlex.quote(c) for c in command))
    result = subprocess.run(command, check=False)
    if result.returncode != 0:
        sys.exit(result.returncode)

    # Uncovered functions per module (terminal, and GitHub job summary if available)
    with open(output / "coverage.json", encoding="utf-8") as stream:
        modules = uncovered_functions.collect(json.load(stream))
    print()
    uncovered_functions.print_summary(modules)
    summary_file = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary_file:
        with open(summary_file, "a", encoding="utf-8") as stream:
            stream.write(uncovered_functions.markdown_summary(modules))
    print(f"\nDetailed report: {output / 'coverage.html'}")
    print("Uncovered functions of a module: python3 "
          f"{source / 'scripts' / 'uncovered_functions.py'} "
          f"{output / 'coverage.json'} -m <Module>")


if __name__ == "__main__":
    main()
