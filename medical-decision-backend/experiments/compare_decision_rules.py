"""
ADDITIONAL ANALYSIS -- ablation of the decision rule, not one of the main
experiments of the manuscript. See README.md in this directory.

Compares the previous decision rule with the current one.

Nothing here alters fusion.py. The previous rule is reimplemented as a
reference so that the two can be placed side by side.

TWO changes are superimposed, and they are separable:

  (a) the circuit      before mu_o = 1 - mu_XOR, nu_o = 1 - nu_XOR (outside U~)
                       now    mu_o = nu_XOR,     nu_o = mu_XOR     (Proposition 2)

  (b) the rule         before q_conf / hesitation / max(1-q_conf, hesitation)
                       now    Xu-Yager score over the mean evidence

The report therefore has THREE columns: previous rule with the previous
degrees (what the program actually did), previous rule with the corrected
degrees (isolating the rule change), and the current rule.

    python3 experiments/compare_decision_rules.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import fusion  # noqa: E402
from qif_xnor import S_P, T_P, xnor_analytic  # noqa: E402


# =====================================================
# Agreement degrees, under both conventions
# =====================================================

def degrees_old(x, y):
    """(1 - mu_XOR, 1 - nu_XOR): the fuzzy negation 1-x, which leaves U~."""
    nu_xnor, mu_xnor = xnor_analytic(x, y)   # (nu_o, mu_o) == (mu_XOR, nu_XOR)
    return 1.0 - mu_xnor, 1.0 - nu_xnor      # returns (mu, nu) in the previous convention


def degrees_new(x, y):
    """Proposition 2."""
    return xnor_analytic(x, y)


# =====================================================
# PREVIOUS rule (reimplementation of the earlier evaluate_rules function)
# =====================================================

def rule_old(x, y, degrees):
    mu_o, nu_o = degrees(x, y)
    pi_o = max(0.0, 1.0 - mu_o - nu_o)

    avg_mu = (x[0] + y[0]) / 2
    avg_nu = (x[1] + y[1]) / 2
    pi_med = (max(0.0, 1 - x[0] - x[1]) + max(0.0, 1 - y[0] - y[1])) / 2

    hat_pi = max(pi_med, pi_o)                 # hesitation
    q_conf = mu_o * (1.0 - nu_o)               # quantum confidence

    s_treat = avg_mu * q_conf * (1.0 - hat_pi)
    s_notreat = avg_nu * q_conf * (1.0 - hat_pi)
    s_exams = max(1.0 - q_conf, hat_pi)

    S = s_treat + s_notreat + s_exams
    if S < 1e-9:
        s_treat, s_notreat, s_exams = 0.33, 0.33, 0.34
    else:
        s_treat, s_notreat, s_exams = s_treat / S, s_notreat / S, s_exams / S

    scores = {"TREAT": s_treat, "DO_NOT_TREAT": s_notreat, "REQUEST_EXAMS": s_exams}
    return {
        "mu_o": mu_o, "nu_o": nu_o, "pi_o": pi_o,
        "sum_munu": mu_o + nu_o,
        "q_conf": q_conf, "hat_pi": hat_pi,
        "scores": scores,
        "decision": max(scores, key=scores.get),
    }


# =====================================================
# CURRENT rule (mirrors fusion.decide for n = 2)
# =====================================================

def rule_new(x, y):
    r = fusion.decide([x, y], exact=True)
    c = r["consensus"]
    return {
        "mu_o": c["mu"], "nu_o": c["nu"], "pi_o": c["pi"],
        "sum_munu": c["mu"] + c["nu"],
        "s_e": r["evidence"]["score"],
        "scores": r["scores"],
        "decision": r["decision"],
    }


CASES = [
    ("agree, treat",        (0.95, 0.02), (0.90, 0.05)),
    ("agree, do not treat", (0.05, 0.90), (0.02, 0.95)),
    ("total conflict",      (0.90, 0.05), (0.05, 0.90)),
    ("UI default",          (0.70, 0.20), (0.65, 0.25)),
    ("moderate agreement",  (0.60, 0.30), (0.55, 0.35)),
    ("high hesitation",     (0.30, 0.30), (0.30, 0.30)),
    ("total ignorance",     (0.00, 0.00), (0.00, 0.00)),
    ("certainty, treat",    (1.00, 0.00), (1.00, 0.00)),
    ("D_I3 identity",       (1.00, 0.00), (0.70, 0.20)),
    ("mild disagreement",   (0.75, 0.15), (0.40, 0.45)),
]


def main():
    print("=" * 108)
    print("AGREEMENT DEGREES".center(108))
    print("=" * 108)
    print(f"{'case':22}{'previous (1-x)':>28}{'':4}{'current (Prop.2)':>26}")
    print(f"{'':22}{'mu_o':>9}{'nu_o':>9}{'soma':>10}{'':4}{'mu_o':>9}{'nu_o':>9}{'pi_o':>8}")
    print("-" * 108)
    for name, x, y in CASES:
        o = rule_old(x, y, degrees_old)
        n = rule_new(x, y)
        flag = " <-- outside U~" if o["sum_munu"] > 1 + 1e-9 else ""
        print(f"{name:22}{o['mu_o']:>9.4f}{o['nu_o']:>9.4f}{o['sum_munu']:>10.4f}"
              f"{'':4}{n['mu_o']:>9.4f}{n['nu_o']:>9.4f}{n['pi_o']:>8.4f}{flag}")

    print()
    print("=" * 108)
    print("INTERMEDIATE VALUES".center(108))
    print("=" * 108)
    print(f"{'case':22}{'PREVIOUS RULE':>34}{'':4}{'CURRENT RULE':>20}")
    print(f"{'':22}{'q_conf':>11}{'hat_pi':>11}{'1-q_conf':>12}{'':4}{'mu_C':>10}{'s(E~)':>10}")
    print("-" * 108)
    for name, x, y in CASES:
        o = rule_old(x, y, degrees_old)
        n = rule_new(x, y)
        print(f"{name:22}{o['q_conf']:>11.4f}{o['hat_pi']:>11.4f}{1-o['q_conf']:>12.4f}"
              f"{'':4}{n['mu_o']:>10.4f}{n['s_e']:>+10.4f}")

    print()
    print("=" * 108)
    print("SCORES AND DECISION".center(108))
    print("=" * 108)
    hdr = f"{'case':22}"
    for tag in ("previous rule + previous degrees", "previous rule + corrected degrees", "current"):
        hdr += f"{tag:>28}"
    print(hdr)
    sub = f"{'':22}"
    for _ in range(3):
        sub += f"{'treat':>7}{'no':>7}{'exam':>7}{'':7}"
    print(sub)
    print("-" * 108)

    changed_by_rule = 0
    for name, x, y in CASES:
        a = rule_old(x, y, degrees_old)
        b = rule_old(x, y, degrees_new)
        c = rule_new(x, y)
        row = f"{name:22}"
        for r in (a, b, c):
            s = r["scores"]
            row += (f"{s['TREAT']:>7.3f}{s['DO_NOT_TREAT']:>7.3f}"
                    f"{s['REQUEST_EXAMS']:>7.3f}{'':7}")
        print(row)
        print(f"{'':22}{a['decision']:>21}{'':7}{b['decision']:>21}{'':7}{c['decision']:>21}")
        if b["decision"] != c["decision"]:
            changed_by_rule += 1

    print("-" * 108)
    print(f"decisions that change because of the RULE (corrected degrees in both): "
          f"{changed_by_rule}/{len(CASES)}")


if __name__ == "__main__":
    main()
