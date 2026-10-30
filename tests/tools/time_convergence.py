"""
Time convergence study: L2 error on phi at the final time vs the time step, for each scheme.

Reads Saves_<scheme>_dt_<dt>/<problem>/time_specialized.csv (4th column: L2 error, last row:
final time), writes <output>.csv (scheme, dt, error) and <output>.png (default output:
time_convergence_output), and prints the fitted slope of each scheme.
With --check, exits with a non-zero code if a slope is lower than the expected order minus
--tolerance (--expected, default "EulerImplicit:1,BDF2:2"; only the listed schemes are checked).
The check is one-sided: it detects a loss of order; a slope above the expected order (error
compensations in the pre-asymptotic regime) is only reported.
"""
import csv
import os
import sys
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_EXPECTED = "EulerImplicit:1,BDF2:2"


def read_results(base_dir, problem):
    results = {}
    for root, _, files in os.walk(base_dir):
        if os.path.basename(root) != problem or "time_specialized.csv" not in files:
            continue
        folder = os.path.basename(os.path.dirname(root))
        if not folder.startswith("Saves_") or "_dt_" not in folder:
            continue
        scheme, dt = folder[len("Saves_"):].rsplit("_dt_", 1)
        df = pd.read_csv(os.path.join(root, "time_specialized.csv"))
        if df.shape[1] < 4:
            continue
        # (dt, error at the last saved time, last saved time)
        results.setdefault(scheme, []).append(
            (float(dt), float(df.iloc[-1, 3]), float(df.iloc[-1, 2])))
    return {k: sorted(v) for k, v in sorted(results.items())}


def main():
    parser = ArgumentParser(description="Time convergence study (EulerImplicit vs BDF2)")
    parser.add_argument("-p", "--problemName", type=str, default="Problem1")
    parser.add_argument("--check", action="store_true", help="check the observed orders")
    parser.add_argument("--tolerance", type=float, default=0.15)
    parser.add_argument("--expected", type=str, default=DEFAULT_EXPECTED,
                        help="expected orders checked by --check, e.g. \"BDF2:2\"")
    parser.add_argument("-o", "--output", type=str, default="time_convergence_output",
                        help="basename of the output files (.csv, .png)")
    args = parser.parse_args()

    expected = {k.strip(): float(v) for k, v in
                (item.split(":") for item in args.expected.split(",") if item.strip())}
    results = read_results("./", args.problemName)
    if not results:
        print("No result found")
        return 1

    status = 0
    # Runs that did not reach the final time (crash, Newton failure...) are excluded: their last
    # error is not at the final time
    t_final = max(t for data in results.values() for _, _, t in data)
    for scheme in list(results):
        complete = []
        for d, e, t in results[scheme]:
            if abs(t - t_final) > 1.e-8 * max(1., abs(t_final)):
                print(f"{scheme}, dt = {d}: INCOMPLETE run (last time {t} < {t_final}), ignored")
                status = 1 if args.check else status
            else:
                complete.append((d, e))
        results[scheme] = complete
    if args.check:
        for scheme in expected:
            if len(results.get(scheme, [])) < 2:
                print(f"{scheme}: FAILED, not enough complete runs to compute an order")
                status = 1
    with open(args.output + ".csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        # Header line (col_compare.py skips the first line of the files it compares)
        writer.writerow(["scheme", "dt", "phi_L2-error"])
        for i, (scheme, data) in enumerate(results.items()):
            if len(data) < 2:
                continue
            dts = np.array([d for d, _ in data])
            errs = np.array([e for _, e in data])
            for d, e in data:
                writer.writerow([scheme, d, e])
            slope, intercept = np.polyfit(np.log(dts), np.log(errs), 1)
            local = np.diff(np.log(errs)) / np.diff(np.log(dts))
            print(f"{scheme}: slope = {slope:.4f}, local orders = {np.round(local, 3)}")
            plt.loglog(dts, errs, "o", color=f"C{i}", label=scheme)
            plt.loglog(dts, np.exp(intercept) * dts**slope, "--", color=f"C{i}",
                       label=f"{scheme} fit: $\\propto \\delta t^{{{slope:.3f}}}$")
            if args.check and scheme in expected:
                if slope < expected[scheme] - args.tolerance:
                    print(f"  -> FAILED: expected order {expected[scheme]}")
                    status = 1
                elif slope > expected[scheme] + args.tolerance:
                    print(f"  -> note: slope above the expected order {expected[scheme]} "
                          "(pre-asymptotic regime)")

    plt.title(f"Time convergence study (MMS, {args.problemName})")
    plt.xlabel(r"Time step $\delta t$")
    plt.ylabel(r"$L_2$ error at $t_f$")
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.savefig(args.output + ".png")
    plt.close()
    return status


if __name__ == "__main__":
    sys.exit(main())
