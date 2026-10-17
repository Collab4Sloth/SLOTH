#!/usr/bin/env python3
"""
Check that the Allen-Cahn interface stays at equilibrium (last time step of
Saves/Problem1/time_specialized.csv):
  - the phase field stays close to the analytical hyperbolic tangent profile
    (normalized L2 error computed by SLOTH);
  - the interfacial energy computed by SLOTH equals sigma (input of main.cpp).
Pure Python (standard library only).
"""
import csv
import sys

RESULTS = "Saves/Problem1/time_specialized.csv"
SIGMA = 6.e-2              # interfacial energy given in main.cpp [J.m-2]
TOL_PROFILE = 1.e-3        # max normalized L2 error with respect to the tanh profile
TOL_SIGMA = 2.e-3          # max relative error on the interfacial energy


def column(row, suffix):
    return float(next(v for k, v in row.items() if k.strip().endswith(suffix)))


def main():
    with open(RESULTS, newline="") as f:
        last = list(csv.DictReader(f))[-1]
    time = column(last, "Time[s]")
    profile_error = column(last, "L2-error normalized[-]")
    sigma = column(last, "Sigma[J.m-3]")
    sigma_error = abs(sigma - SIGMA) / SIGMA

    ok_profile = profile_error < TOL_PROFILE
    ok_sigma = sigma_error < TOL_SIGMA
    print(f"t = {time} s")
    print(f"normalized L2 error / tanh profile = {profile_error:.3e} (< {TOL_PROFILE}) "
          f"{'OK' if ok_profile else 'FAILED'}")
    print(f"interfacial energy = {sigma:.5e} (sigma = {SIGMA}, relative error {sigma_error:.2e} "
          f"< {TOL_SIGMA}) {'OK' if ok_sigma else 'FAILED'}")
    return 0 if (ok_profile and ok_sigma) else 1


if __name__ == "__main__":
    sys.exit(main())
