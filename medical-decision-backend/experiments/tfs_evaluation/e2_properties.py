"""
E2 -- Numerical validation of the properties of boxplus_IP.

NUMERICAL verification over grids of U~; it does not replace the analytical
proof.
"""

from typing import Dict, List

from _common import (GRID_STEP_CLOSURE, GRID_STEP_PROPERTIES, SEEDS, TOL_EXACT,
                     TOL_SUM, banner, sample_opinions, valid_grid, write_csv)

from fusion import aggregate
from qif_xnor import IFV_ONE, IFV_ZERO, xnor_analytic

ROWS: List[Dict] = []


def record(prop, domain, n_inputs, n_pairs, violations, max_err, tol, verdict):
    ROWS.append({
        "property": prop, "domain": domain, "n_inputs": n_inputs,
        "n_pairs": n_pairs, "violations": violations,
        "max_abs_error": f"{max_err:.3e}", "tolerance": f"{tol:.1e}",
        "verdict": verdict,
    })


def main() -> List[Dict]:
    banner("E2  XNOR PROPERTIES (numerical verification)")
    gp = valid_grid(GRID_STEP_PROPERTIES)
    gc = valid_grid(GRID_STEP_CLOSURE)
    print(f"  property grid: step {GRID_STEP_PROPERTIES}"
          f" -> {len(gp)} valid points, {len(gp)**2} ordered pairs")
    print(f"  closure grid : step {GRID_STEP_CLOSURE}"
          f" -> {len(gc)} valid points, {len(gc)**2} ordered pairs")

    checks = [
        (IFV_ZERO, IFV_ZERO, (1.0, 0.0)),
        (IFV_ONE, IFV_ONE, (1.0, 0.0)),
        (IFV_ZERO, IFV_ONE, (0.0, 1.0)),
        (IFV_ONE, IFV_ZERO, (0.0, 1.0)),
    ]
    v = 0
    err = 0.0
    for a, b, exp in checks:
        got = xnor_analytic(a, b)
        e = max(abs(got[0] - exp[0]), abs(got[1] - exp[1]))
        err = max(err, e)
        if e > TOL_EXACT:
            v += 1
    record("D_I1 boundary", "4 extreme combinations", 2, len(checks), v, err,
           TOL_EXACT, "pass" if v == 0 else "FAIL")
    print(f"  D_I1  {len(checks)} combinations  violations {v}  max error {err:.3e}")

    v = 0
    err = 0.0
    for x in gp:
        for y in gp:
            a, b = xnor_analytic(x, y), xnor_analytic(y, x)
            e = max(abs(a[0] - b[0]), abs(a[1] - b[1]))
            err = max(err, e)
            if e > TOL_EXACT:
                v += 1
    record("D_I2 symmetry", f"grid step {GRID_STEP_PROPERTIES}", len(gp),
           len(gp) ** 2, v, err, TOL_EXACT, "pass" if v == 0 else "FAIL")
    print(f"  D_I2  {len(gp)**2} pairs  violations {v}  max error {err:.3e}")

    v = 0
    err = 0.0
    for x in gc:
        got = xnor_analytic(IFV_ONE, x)
        e = max(abs(got[0] - x[0]), abs(got[1] - x[1]))
        err = max(err, e)
        if e > TOL_EXACT:
            v += 1
    record("D_I3 identity", f"grid step {GRID_STEP_CLOSURE}", len(gc), len(gc),
           v, err, TOL_EXACT, "pass" if v == 0 else "FAIL")
    print(f"  D_I3  {len(gc)} values  violations {v}  max error {err:.3e}")

    v_neg = v_sum = 0
    worst = 0.0
    for x in gc:
        for y in gc:
            mu, nu = xnor_analytic(x, y)
            if mu < -TOL_EXACT or nu < -TOL_EXACT:
                v_neg += 1
            excess = mu + nu - 1.0
            worst = max(worst, excess)
            if excess > TOL_SUM:
                v_sum += 1
    record("closure mu >= 0 and nu >= 0", f"grid step {GRID_STEP_CLOSURE}", len(gc),
           len(gc) ** 2, v_neg, 0.0, TOL_EXACT, "pass" if v_neg == 0 else "FAIL")
    record("closure mu+nu <= 1", f"grid step {GRID_STEP_CLOSURE}", len(gc),
           len(gc) ** 2, v_sum, max(0.0, worst), TOL_SUM,
           "pass" if v_sum == 0 else "FAIL")
    print(f"  closure  {len(gc)**2} pairs  non-negativity {v_neg}"
          f"  sum<=1 {v_sum}  largest excess {worst:.3e}")

    for n in (2, 3, 4, 5, 8, 10, 12):
        v = 0
        worst_n = 0.0
        for rep in range(200):
            ops = sample_opinions(n, SEEDS["opinion_sampling"] + 1000 * n + rep)
            c, _, _ = aggregate(ops, exact=True)
            worst_n = max(worst_n, c[0] + c[1] - 1.0)
            if c[0] < -TOL_EXACT or c[1] < -TOL_EXACT or c[0] + c[1] > 1 + TOL_SUM:
                v += 1
        record(f"multi-source closure n={n}", "200 random valid opinion sets",
               n, n * (n - 1) // 2, v, max(0.0, worst_n), TOL_SUM,
               "pass" if v == 0 else "FAIL")
        print(f"  multi-source n={n:2}  200 sets  violations {v}"
              f"  largest excess {worst_n:.3e}")

    path = write_csv("raw_property_validation.csv", list(ROWS[0]), ROWS)
    print(f"\n  -> {path.rsplit('/', 1)[-1]}")
    print(f"  total rows: {len(ROWS)};  "
          f"all 'pass': {all(r['verdict'] == 'pass' for r in ROWS)}")
    return ROWS


if __name__ == "__main__":
    main()
