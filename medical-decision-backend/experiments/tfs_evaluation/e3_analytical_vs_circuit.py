"""
E3 -- Agreement between the closed analytical form and the sampled circuit.

For each representative case and each shot budget, compares boxplus_IP
evaluated in closed form with the estimate obtained by measurement.
"""

import math
from typing import Dict, List

from _common import (PRODUCTION_DEFAULT_SHOTS, REPRESENTATIVE_CASES, SEEDS,
                     SHOT_BUDGETS, banner, write_csv)

from qif_xnor import run_xnor_joint, xnor_analytic

# The production default is included in the representative-case table.
ALL_BUDGETS = tuple(sorted(set(SHOT_BUDGETS) | {PRODUCTION_DEFAULT_SHOTS}))

# A single draw per (case, budget) gives a MAE too noisy to support any
# claim about the decay rate. The summary uses R independent repetitions
# per cell.
R_REPS = 30


def case_seed(case_index: int, budget_index: int) -> int:
    """Deterministic seed per (case, budget). Documented in the README."""
    return SEEDS["analytical_vs_circuit"] + 100 * case_index + budget_index


def main():
    banner("E3  ANALYTICAL vs QUANTUM-CIRCUIT AGREEMENT")
    print(f"  shot budgets: {ALL_BUDGETS}  (production default = "
          f"{PRODUCTION_DEFAULT_SHOTS})")

    rows: List[Dict] = []
    for ci, (regime, x, y) in enumerate(REPRESENTATIVE_CASES):
        mu_a, nu_a = xnor_analytic(x, y)
        pi_a = 1.0 - mu_a - nu_a
        for bi, shots in enumerate(ALL_BUDGETS):
            seed = case_seed(ci, bi)
            j = run_xnor_joint(x, y, shots=shots, seed=seed)
            mu_c, nu_c = j["mu"], j["nu"]
            pi_c = 1.0 - mu_c - nu_c
            rows.append({
                "regime": regime,
                "mu_x": x[0], "nu_x": x[1], "mu_y": y[0], "nu_y": y[1],
                "shots": shots, "seed": seed,
                "mu_analytical": f"{mu_a:.6f}", "mu_circuit": f"{mu_c:.6f}",
                "mu_abs_error": f"{abs(mu_c - mu_a):.6f}",
                "nu_analytical": f"{nu_a:.6f}", "nu_circuit": f"{nu_c:.6f}",
                "nu_abs_error": f"{abs(nu_c - nu_a):.6f}",
                "pi_analytical": f"{pi_a:.6f}", "pi_circuit": f"{pi_c:.6f}",
                "pi_abs_error": f"{abs(pi_c - pi_a):.6f}",
            })

    write_csv("analytical_vs_circuit.csv", list(rows[0]), rows)

    print(f"\n  cases at the production default ({PRODUCTION_DEFAULT_SHOTS} shots):")
    print(f"  {'regime':32}{'mu_an':>9}{'mu_ci':>9}{'|e|':>9}"
          f"{'nu_an':>9}{'nu_ci':>9}{'|e|':>9}")
    for r in rows:
        if r["shots"] == PRODUCTION_DEFAULT_SHOTS:
            print(f"  {r['regime']:32}{r['mu_analytical']:>9}{r['mu_circuit']:>9}"
                  f"{r['mu_abs_error']:>9}{r['nu_analytical']:>9}"
                  f"{r['nu_circuit']:>9}{r['nu_abs_error']:>9}")

    # ---- summary per budget, with R repetitions per cell ----
    summary: List[Dict] = []
    print(f"\n  summary over {len(REPRESENTATIVE_CASES)} cases x {R_REPS}"
          f" independent repetitions ({len(REPRESENTATIVE_CASES) * R_REPS}"
          f" estimates per budget):")
    print(f"  {'shots':>8}{'MAE mu':>11}{'MAE nu':>11}{'MAE pi':>11}"
          f"{'RMSE mu':>11}{'RMSE nu':>11}{'max |e|':>11}{'sqrt(N)*MAE':>13}")
    for bi, shots in enumerate(ALL_BUDGETS):
        emu, enu, epi = [], [], []
        for ci, (regime, x, y) in enumerate(REPRESENTATIVE_CASES):
            mu_a, nu_a = xnor_analytic(x, y)
            pi_a = 1.0 - mu_a - nu_a
            for rep in range(R_REPS):
                sd = case_seed(ci, bi) + 10_000 * (rep + 1)
                j = run_xnor_joint(x, y, shots=shots, seed=sd)
                emu.append(abs(j["mu"] - mu_a))
                enu.append(abs(j["nu"] - nu_a))
                epi.append(abs((1.0 - j["mu"] - j["nu"]) - pi_a))
        mae_mu = sum(emu) / len(emu)
        mae_nu = sum(enu) / len(enu)
        mae_pi = sum(epi) / len(epi)
        rmse_mu = math.sqrt(sum(e * e for e in emu) / len(emu))
        rmse_nu = math.sqrt(sum(e * e for e in enu) / len(enu))
        mx = max(emu + enu + epi)
        scaled = math.sqrt(shots) * mae_mu       # deve ficar ~constante
        summary.append({
            "shots": shots, "n_cases": len(REPRESENTATIVE_CASES),
            "n_repetitions": R_REPS, "n_estimates": len(emu),
            "mae_mu": f"{mae_mu:.6f}", "mae_nu": f"{mae_nu:.6f}",
            "mae_pi": f"{mae_pi:.6f}",
            "rmse_mu": f"{rmse_mu:.6f}", "rmse_nu": f"{rmse_nu:.6f}",
            "max_abs_error": f"{mx:.6f}",
            "sqrt_shots_times_mae_mu": f"{scaled:.4f}",
        })
        print(f"  {shots:>8}{mae_mu:>11.6f}{mae_nu:>11.6f}{mae_pi:>11.6f}"
              f"{rmse_mu:>11.6f}{rmse_nu:>11.6f}{mx:>11.6f}{scaled:>13.4f}")

    write_csv("shot_summary.csv", list(summary[0]), summary)
    print("\n  -> analytical_vs_circuit.csv, shot_summary.csv")
    print("  The sqrt(N)*MAE column would be constant if the error decayed as"
          " 1/sqrt(N).")
    return rows, summary


if __name__ == "__main__":
    main()
