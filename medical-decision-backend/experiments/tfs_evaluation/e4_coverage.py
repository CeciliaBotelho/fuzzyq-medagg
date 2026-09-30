"""
E4 -- Finite-shot analysis: empirical coverage of the 95% intervals.

For each configuration and each shot budget, runs R independent repetitions
and measures how often the reported interval contains the exact analytical
value.

Saves the per-repetition result (raw) and the aggregated summary.
"""

import statistics
import time
from typing import Dict, List

from _common import SEEDS, SHOT_BUDGETS, banner, write_csv

from fusion import aggregate, decide

R_REPS = 200

# Multi-source configurations, with the actual tuples used.
CONFIGS = [
    ("n=2", [(0.90, 0.05), (0.85, 0.10)]),
    ("n=3", [(0.90, 0.05), (0.85, 0.10), (0.20, 0.70)]),
    ("n=4", [(0.90, 0.05), (0.85, 0.10), (0.20, 0.70), (0.55, 0.35)]),
]


def main():
    banner("E4  FINITE-SHOT ANALYSIS: empirical coverage of the 95% intervals")
    print(f"  {R_REPS} independent repetitions per cell;"
          f" shot budgets {SHOT_BUDGETS}")
    print(f"  seed of repetition r = {SEEDS['coverage']}"
          f" + 1000*config_index + r")

    raw: List[Dict] = []
    summary: List[Dict] = []

    for ci, (label, ops) in enumerate(CONFIGS):
        truth, _, _ = aggregate(ops, exact=True)
        mu_t, nu_t = truth
        pi_t = 1.0 - mu_t - nu_t
        n_pairs = len(ops) * (len(ops) - 1) // 2
        print(f"\n  {label}: {n_pairs} pairs, analytical truth"
              f" mu={mu_t:.6f} nu={nu_t:.6f} pi={pi_t:.6f}")

        for shots in SHOT_BUDGETS:
            t0 = time.perf_counter()
            hits = {"mu": 0, "nu": 0, "pi": 0}
            hw = {"mu": [], "nu": [], "pi": []}

            for rep in range(R_REPS):
                seed = SEEDS["coverage"] + 1000 * ci + rep
                r = decide(ops, shots=shots, seed=seed)
                c, u = r["consensus"], r["uncertainty"]
                cov = {
                    "mu": abs(c["mu"] - mu_t) <= u["mu_ci95"],
                    "nu": abs(c["nu"] - nu_t) <= u["nu_ci95"],
                    "pi": abs(c["pi"] - pi_t) <= u["pi_ci95"],
                }
                for k in hits:
                    hits[k] += int(cov[k])
                hw["mu"].append(u["mu_ci95"])
                hw["nu"].append(u["nu_ci95"])
                hw["pi"].append(u["pi_ci95"])

                raw.append({
                    "config": label, "n_opinions": len(ops), "n_pairs": n_pairs,
                    "shots": shots, "seed": seed, "repetition": rep,
                    "mu_truth": f"{mu_t:.8f}", "nu_truth": f"{nu_t:.8f}",
                    "pi_truth": f"{pi_t:.8f}",
                    "mu_hat": f"{c['mu']:.8f}", "nu_hat": f"{c['nu']:.8f}",
                    "pi_hat": f"{c['pi']:.8f}",
                    "mu_ci95": f"{u['mu_ci95']:.8f}",
                    "nu_ci95": f"{u['nu_ci95']:.8f}",
                    "pi_ci95": f"{u['pi_ci95']:.8f}",
                    "mu_covered": int(cov["mu"]), "nu_covered": int(cov["nu"]),
                    "pi_covered": int(cov["pi"]),
                })

            dt = time.perf_counter() - t0
            row = {
                "config": label, "n_opinions": len(ops), "n_pairs": n_pairs,
                "shots": shots, "repetitions": R_REPS,
                "coverage_mu": f"{hits['mu'] / R_REPS:.3f}",
                "coverage_nu": f"{hits['nu'] / R_REPS:.3f}",
                "coverage_pi": f"{hits['pi'] / R_REPS:.3f}",
                "mean_halfwidth_mu": f"{statistics.fmean(hw['mu']):.6f}",
                "mean_halfwidth_nu": f"{statistics.fmean(hw['nu']):.6f}",
                "mean_halfwidth_pi": f"{statistics.fmean(hw['pi']):.6f}",
                "median_halfwidth_mu": f"{statistics.median(hw['mu']):.6f}",
                "median_halfwidth_nu": f"{statistics.median(hw['nu']):.6f}",
                "median_halfwidth_pi": f"{statistics.median(hw['pi']):.6f}",
                "wall_seconds": f"{dt:.1f}",
            }
            summary.append(row)
            print(f"    shots {shots:>6}: coverage mu={row['coverage_mu']}"
                  f" nu={row['coverage_nu']} pi={row['coverage_pi']}"
                  f"  half-width mu={row['mean_halfwidth_mu']}"
                  f"  ({dt:.0f}s)")

    write_csv("finite_shot_coverage.csv", list(summary[0]), summary)
    write_csv("finite_shot_coverage_raw.csv", list(raw[0]), raw)

    allc = [float(r[f"coverage_{k}"]) for r in summary for k in ("mu", "nu", "pi")]
    print(f"\n  {len(summary)} cells, {len(raw)} repetitions in total")
    print(f"  minimum coverage {min(allc):.3f}  maximum {max(allc):.3f}"
          f"  mean {statistics.fmean(allc):.4f}   (nominal 0.95)")
    print("  -> finite_shot_coverage.csv, finite_shot_coverage_raw.csv")
    return summary


if __name__ == "__main__":
    main()
