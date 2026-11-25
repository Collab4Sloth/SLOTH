#!/usr/bin/env python3
"""
Outputs of PFHub benchmark 3 from the result folders of DendriticGrowth2Dtest1 (option -o):
tip position along x (largest x of the isovalue phi = 0), solid fraction (from phi_integral,
int (phi + 1) / 2 divided by the area of the domain), free energy (Density) as functions of time.
Writes <folder>/benchmark3.csv for each folder and, with matplotlib, benchmark3.png comparing the
folders with the PFHub submissions of pfhub_bm3_references.csv (MOOSE, AMPE, MEMPHIS; data of
https://github.com/usnistgov/pfhub, _data/simulations/*/meta.yaml).

  python3 plot.py Saves_euler_dt0.01 Saves_bdf2_dt0.04 [-L 960] [--no-pfhub]
"""
import csv
import sys
from argparse import ArgumentParser
from pathlib import Path


def read_columns(path):
    """Rows of time_specialized.csv: the last run if several runs were appended"""
    with open(path) as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = []
        for row in reader:
            if rows and int(row[0]) <= int(rows[-1][0]):
                rows = []  # new run appended to the file
            rows.append(row)
    return header, rows


def tip_positions(path):
    """Time -> largest x and y of the isovalue (last run of the file)"""
    tips = {}
    last = -1
    with open(path) as f:
        f.readline()
        for line in f:
            v = line.rstrip("\n").rstrip(",").split(",")
            it = int(v[0])
            if it <= last:
                tips = {}
            last = it
            xs = [float(x) for x in v[3::2] if x]
            ys = [float(y) for y in v[4::2] if y]
            if xs:
                tips[float(v[2])] = (max(xs), max(ys))
    return tips


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("folders", nargs="+")
    parser.add_argument("-L", "--length", type=float, default=960.,
                        help="size of the square domain (default 960)")
    parser.add_argument("--no-pfhub", action="store_true",
                        help="do not plot the PFHub submissions")
    args = parser.parse_args()
    area = args.length * args.length
    curves = {}
    for folder in args.folders:
        pf = Path(folder) / "PhaseField"
        header, rows = read_columns(pf / "time_specialized.csv")
        i_t, i_int = header.index("Time[s]"), header.index("phi_integral[-]")
        i_f = header.index("Density[J.m-3]")
        tips = tip_positions(pf / "phi_iso_computation.csv")
        out = []
        for row in rows:
            t = float(row[i_t])
            solid = 0.5 * (float(row[i_int]) + area) / area  # int (phi + 1) / 2 / area
            tip_x, tip_y = tips.get(t, (float("nan"), float("nan")))
            out.append((t, solid, float(row[i_f]), tip_x, tip_y))
        with open(pf / "benchmark3.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["time", "solid_fraction", "free_energy", "tip_x", "tip_y"])
            writer.writerows(out)
        curves[folder] = out
        t, solid, energy, tip_x, _ = out[-1]
        print(f"{folder}: t = {t:g}, solid fraction = {solid:.6g}, free energy = {energy:.6g}, "
              f"tip x = {tip_x:.4g}")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return 0
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    references = Path(__file__).resolve().parent / "pfhub_bm3_references.csv"
    if not args.no_pfhub and references.exists():
        series = {}
        with open(references) as f:
            for row in csv.DictReader(f):
                series.setdefault((row["simulation"], row["quantity"]), []).append(
                    (float(row["time"]), float(row["value"])))
        labelled = False
        for (_, quantity), points in series.items():
            k = {"tip_position": 0, "solid_fraction": 1, "free_energy": 2}[quantity]
            points.sort()
            label = None
            if k == 0 and not labelled:
                label, labelled = "PFHub submissions (MOOSE, AMPE, MEMPHIS)", True
            axes[k].plot([p[0] for p in points], [p[1] for p in points], color="#9a9993",
                         linestyle="--", linewidth=1, label=label)
    for folder, out in curves.items():
        t = [o[0] for o in out]
        for ax, k in zip(axes, (3, 1, 2)):
            ax.plot(t, [o[k] for o in out], label=folder)
    for ax, title in zip(axes, ("tip position (x)", "solid fraction", "free energy")):
        ax.set_xlabel("time")
        ax.set_title(title)
        ax.grid(True, color="#e6e5e0")
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig("benchmark3.png", dpi=120)
    return 0


if __name__ == "__main__":
    sys.exit(main())
