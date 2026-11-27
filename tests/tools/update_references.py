#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check and update the reference files of the comparison tests (col_compare.py) of SLOTH.

For each comparison test of a build directory (tests created by create_col_comparison and
create_col_comparison_convergence), the results of the last run are compared with the reference,
column by column, with the criterion of the test. The script reports the columns which differ and,
with --update, copies the results on the reference of the source tree (tests/<...>/ref/<file>)
and on its copy in the build directory (ref_<test>/<file>), so that the tests pass without
configuring again.

Run the tests first (ctest), then, from anywhere:

  # report of all the comparison tests (nothing is modified)
  python3 tests/tools/update_references.py -b build

  # tests which failed during the last ctest run, whose name matches a regular expression
  python3 tests/tools/update_references.py -b build --failed -R AllenCahn

  # update the references which differ in the column Density[J.m-3] only
  python3 tests/tools/update_references.py -b build --update --columns "Density[J.m-3]"

Only the Python standard library is used (and col_compare.py, next to this script).
"""

import os
import re
import shutil
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from col_compare import calculate_error, read_csv, _to_float  # noqa: E402


class Comparison:
    """Comparison test: files, columns and criterion of its col_compare.py command"""

    def __init__(self, name, directory, results, reference, columns, error, criterion):
        self.name = name
        self.directory = directory  # build directory of the test (working directory)
        self.results = results
        self.reference = reference  # copy of the reference in the build directory
        self.columns = columns
        self.error = error
        self.criterion = criterion


def quoted_arguments(text):
    """Arguments of a command written in a CTestTestfile.cmake: quoted strings"""
    return re.findall(r'"((?:[^"\\]|\\.)*)"', text)


def comparison_of(name, directory, arguments):
    """Comparison of a test from its arguments, None if it is not a col_compare.py command"""
    # The convergence comparisons run col_compare.py inside a "bash -c" command
    words = []
    for argument in arguments:
        words.extend(argument.replace('\\"', '"').replace('"', ' ').split())
    if "col_compare.py" not in words:
        return None
    try:
        i = words.index("-f")
        results, reference = words[i + 1], words[i + 2]
        error = words[words.index("-e") + 1]
        criterion = float(words[words.index("-k") + 1])
        j = words.index("-c") + 1
        columns = []
        while j < len(words) and re.fullmatch(r"-?\d+", words[j]):
            columns.append(int(words[j]))
            j += 1
    except (ValueError, IndexError):
        return None
    results = Path(results) if os.path.isabs(results) else directory / results
    return Comparison(name, directory, results, Path(reference), columns, error, criterion)


def test_files(directory):
    """CTestTestfile.cmake of the directory and of its subdirectories listed by subdirs()"""
    testfile = directory / "CTestTestfile.cmake"
    if not testfile.exists():
        return []
    files = [testfile]
    for subdirectory in re.findall(r'^subdirs\("([^"]*)"\)', testfile.read_text(), re.MULTILINE):
        files.extend(test_files(directory / subdirectory))
    return files


def comparisons_of(build):
    """All the comparison tests of the build directory"""
    found = []
    for testfile in test_files(build):
        text = testfile.read_text()
        for match in re.finditer(r"^add_test\(\s*(?:\[=\[)?([\w.\-]+)(?:\]=\])?\s+(.*)\)\s*$",
                                 text, re.MULTILINE):
            comparison = comparison_of(match.group(1), testfile.parent,
                                       quoted_arguments(match.group(2)))
            if comparison is not None:
                found.append(comparison)
    return found


def failed_tests(build):
    """Names of the tests which failed during the last ctest run"""
    log = build / "Testing" / "Temporary" / "LastTestsFailed.log"
    if not log.exists():
        return set()
    return {line.split(":", 1)[1].strip() for line in log.read_text().splitlines() if ":" in line}


def source_reference(comparison, build, source):
    """Reference of the source tree: <source>/<test directory>/ref/<file>"""
    directory = comparison.reference.parent
    if not directory.name.startswith("ref_"):
        # reference in a subdirectory of ref_<test>: ref_<test>/<subdirectory>/<file>
        parts = comparison.reference.parts
        index = max(i for i, part in enumerate(parts) if part.startswith("ref_"))
        directory = Path(*parts[:index + 1])
        relative = Path(*parts[index + 1:])
    else:
        relative = Path(comparison.reference.name)
    test_directory = directory.parent.relative_to(build)
    return source / test_directory / "ref" / relative


def differences(comparison):
    """
    (status, columns): status 'ok', 'diff', 'missing results', 'missing reference' or
    'rows'; columns: largest error of each column which differs (as col_compare.py)
    """
    if not comparison.results.exists():
        return "missing results", {}
    if not comparison.reference.exists():
        return "missing reference", {}
    headers, rows = read_csv(comparison.results)
    _, reference_rows = read_csv(comparison.reference)
    columns = comparison.columns
    if not columns or columns[0] < 0:
        columns = list(range(len(headers)))
    errors = {}
    # The first listed column (iteration/time) is not compared, as in col_compare.py
    for column in columns[1:]:
        for row, reference_row in zip(rows, reference_rows):
            if column >= len(row) or column >= len(reference_row):
                continue
            try:
                error = calculate_error(_to_float(row[column]), _to_float(reference_row[column]),
                                        comparison.error)
            except ValueError:
                continue  # text column
            if error > comparison.criterion:
                name = headers[column]
                errors[name] = max(errors.get(name, 0.), error)
    if len(rows) != len(reference_rows):
        return "rows", errors
    return ("diff" if errors else "ok"), errors


def main():
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-b", "--build", required=True, help="build directory (tests already run)")
    parser.add_argument("-s", "--source", help="source directory (default: from CMakeCache.txt)")
    parser.add_argument("-R", "--regex", help="only the tests whose name matches this expression")
    parser.add_argument("--failed", action="store_true",
                        help="only the tests which failed during the last ctest run")
    parser.add_argument("--update", action="store_true",
                        help="copy the results on the references of the tests which differ")
    parser.add_argument("--columns", nargs="+",
                        help="with --update: only the tests whose differences are in these columns")
    parser.add_argument("--force", action="store_true",
                        help="with --update: also the tests whose number of rows differs")
    parser.add_argument("--prefix-map", metavar="OLD=NEW",
                        help="replace the prefix OLD of the paths written by CMake by NEW (build "
                             "configured in another directory)")
    args = parser.parse_args()

    build = Path(args.build).resolve()
    old, new = args.prefix_map.split("=", 1) if args.prefix_map else (None, None)
    if args.source:
        source = Path(args.source).resolve()
    else:
        cache = (build / "CMakeCache.txt").read_text()
        source = Path(re.search(r"^CMAKE_HOME_DIRECTORY:INTERNAL=(.*)$", cache,
                                re.MULTILINE).group(1))
        if old is not None:
            source = Path(str(source).replace(old, new, 1))

    comparisons = comparisons_of(build)
    if old is not None:
        for c in comparisons:
            c.reference = Path(str(c.reference).replace(old, new, 1))
            c.results = Path(str(c.results).replace(old, new, 1))
    if args.regex:
        comparisons = [c for c in comparisons if re.search(args.regex, c.name)]
    if args.failed:
        failed = failed_tests(build)
        comparisons = [c for c in comparisons if c.name in failed]

    counts = {}
    updated = 0
    for comparison in comparisons:
        status, errors = differences(comparison)
        counts[status] = counts.get(status, 0) + 1
        if status == "ok":
            continue
        details = ", ".join(f"{name}: {error:.1e}" for name, error in errors.items())
        print(f"{comparison.name:55s} {status:17s} {details}")
        if not args.update or status in ("missing results", "missing reference"):
            continue
        if status == "rows" and not args.force:
            print(f"{'':55s} not updated (number of rows differs, see --force)")
            continue
        if args.columns and not set(errors) <= set(args.columns):
            print(f"{'':55s} not updated (other columns than {args.columns})")
            continue
        target = source_reference(comparison, build, source)
        if not target.exists():
            print(f"{'':55s} not updated (reference not found in the sources: {target})")
            continue
        shutil.copyfile(comparison.results, target)
        shutil.copyfile(comparison.results, comparison.reference)
        updated += 1
        print(f"{'':55s} updated: {target.relative_to(source)}")

    summary = ", ".join(f"{count} {status}" for status, count in sorted(counts.items()))
    print(f"\n{len(comparisons)} comparison tests: {summary or 'none'}"
          + (f"; {updated} references updated" if args.update else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
