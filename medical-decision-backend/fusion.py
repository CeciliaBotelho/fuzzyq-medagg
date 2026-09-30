"""
Evidence Fusion Core
====================

Aggregates n intuitionistic fuzzy opinions produced by AI models or by
experts, possibly distributed across institutions, and maps the result onto a
clinical action.

Two deliberately separated stages.

1. AGREEMENT -- given by the intuitionistic fuzzy XNOR.

   boxplus_I is a binary EQUIVALENCE operator: both of its operands must be
   opinions. Applying it as a sequential fold breaks that condition from
   n = 3 on, because boxplus_I(x~_1, x~_2) is already an agreement value and
   not an opinion; the fold then stops being direction-agnostic, which is the
   property characterising an equivalence. boxplus_I is also non-associative,
   so no fold is available in the first place.

   Agreement is therefore formed from the P = C(n,2) pairs:

       (mu_ij, nu_ij) = boxplus_I(x~_i, x~_j),   i < j

       C~ = ( SUM w_ij mu_ij ,  SUM w_ij nu_ij ),   SUM w_ij = 1

   Since C~ is a convex combination of points of U~ and U~ is convex, C~ lies
   in U~. And since the weights depend only on the unordered PAIR, C~ does
   not depend on the order in which the sources are presented.

   The method requires n >= 2 and uses a single pipeline: there is no
   separate branch for n = 2. In that case P = C(2,2) = 1, the only pair is
   (1,2), w_12 = 1, and the aggregate reduces exactly to
   boxplus_I(x~_1, x~_2).

2. DIRECTION -- the XNOR is direction-agnostic by construction: two sources
   agreeing on "treat" and two agreeing on "do not treat" produce the same
   mu_C. Direction comes from the score s(x~) = x1 - x2, the same function
   used in the Xu-Yager total order, evaluated on the mean evidence E~.

   The mapping from s to the scores below is NOT derived from the Xu-Yager
   order: it is the decision layer defined by FuzzyQ-MedAgg. The score only
   supplies the direction of the evidence; the rule that converts it into a
   clinical action is a choice of the framework, external to the theory of
   the XNOR.

Weighting schemes
-----------------
uniform       w_ij = 1/P.
weighted      w_ij proportional to w_i * w_j, with per-source weights.
hierarchical  omega_i = 1/(K * n_k) for a source of institution k, among K
              institutions with n_k sources each. The SOURCE weights are
              balanced per institution -- each sums to 1/K -- before the pair
              weights are induced. This does NOT mean equal weight in the
              final aggregate: the intra-institutional mass is
              (n_k-1)/(2 K^2 n_k) and grows with n_k, so a larger site
              retains a somewhat larger share. The domination effect is
              strongly attenuated, not eliminated. Reduces to the weighted
              case.

Uncertainty
-----------
Each mu_ij is a binomial proportion estimated from `shots` samples, and the
pairs use independent random streams. Hence

    Var(mu_C) = SUM w_ij^2 * mu_ij (1 - mu_ij) / shots

and analogously for nu_C. Since mu and nu are read from the SAME measurement
they are correlated, and

    Var(pi_C) = Var(mu_C) + Var(nu_C) + 2 Cov(mu_C, nu_C),
    Cov_ij    = (p11_ij - mu_ij nu_ij) / shots.
"""

import itertools
import math
from typing import Dict, List, Optional, Sequence, Tuple

from qif_xnor import IFV, S_P, T_P, derive_seeds, provenance, run_xnor_joint, xnor_analytic

__all__ = [
    "InvalidOpinion",
    "validate_opinion",
    "mean_ifv",
    "score",
    "accuracy",
    "T_I_P",
    "pair_weights",
    "aggregate",
    "decide",
]

DECISION_TREAT = "TREAT"
DECISION_DO_NOT_TREAT = "DO_NOT_TREAT"
DECISION_REQUEST_EXAMS = "REQUEST_EXAMS"

COMBINE_MODES = ("uniform", "weighted", "hierarchical", "tnorm")

_TOL = 1e-9
_Z95 = 1.959963984540054


class InvalidOpinion(ValueError):
    """Input outside U~, or inconsistent with the given weights/groups."""


# =====================================================
# 1) VALIDATION IN U~ (Eq. 1)
# =====================================================

def validate_opinion(mu: float, nu: float, label: str = "opinion") -> IFV:
    if not (0.0 <= mu <= 1.0):
        raise InvalidOpinion(f"{label}: mu={mu} is outside [0,1]")
    if not (0.0 <= nu <= 1.0):
        raise InvalidOpinion(f"{label}: nu={nu} is outside [0,1]")
    if mu + nu > 1.0 + _TOL:
        raise InvalidOpinion(
            f"{label}: mu+nu={mu + nu:.4f} > 1 violates the intuitionistic condition"
        )
    return float(mu), float(nu)


# =====================================================
# 2) XU-YAGER ORDER (Section III-A)
# =====================================================

def score(x: IFV) -> float:
    """s(x~) = x1 - x2, in [-1, 1]."""
    return x[0] - x[1]


def accuracy(x: IFV) -> float:
    """h(x~) = x1 + x2, in [0, 1]."""
    return x[0] + x[1]


def mean_ifv(opinions: Sequence[IFV], weights: Optional[Sequence[float]] = None) -> IFV:
    """Mean, optionally weighted. U~ is convex, so the result stays in U~."""
    if weights is None:
        n = len(opinions)
        return (sum(o[0] for o in opinions) / n, sum(o[1] for o in opinions) / n)
    total = sum(weights)
    if total <= 0:
        raise InvalidOpinion("weights must sum to a positive value")
    return (
        sum(w * o[0] for w, o in zip(weights, opinions)) / total,
        sum(w * o[1] for w, o in zip(weights, opinions)) / total,
    )


def T_I_P(x: IFV, y: IFV) -> IFV:
    """Representable intersection T_I under T_P/S_P (Eq. repT). Closed in U~."""
    return T_P(x[0], y[0]), S_P(x[1], y[1])


# =====================================================
# 3) PAIR WEIGHTS
# =====================================================

def pair_weights(
    n: int,
    combine: str = "uniform",
    weights: Optional[Sequence[float]] = None,
    groups: Optional[Sequence[object]] = None,
) -> Dict[Tuple[int, int], float]:
    """
    Returns {(i,j): w_ij} over the pairs i < j, with SUM w_ij = 1.

    The weights depend only on the unordered pair, which makes the
    aggregation independent of the order in which the sources are presented.
    """
    pairs = list(itertools.combinations(range(n), 2))
    if not pairs:
        return {}

    if combine == "uniform":
        raw = {p: 1.0 for p in pairs}

    elif combine == "weighted":
        if weights is None:
            raise InvalidOpinion("combine='weighted' requires per-source weights")
        if len(weights) != n:
            raise InvalidOpinion(f"expected {n} weights, got {len(weights)}")
        if any(w < 0 for w in weights):
            raise InvalidOpinion("weights must be non-negative")
        raw = {(i, j): weights[i] * weights[j] for i, j in pairs}

    elif combine == "hierarchical":
        if groups is None:
            raise InvalidOpinion("combine='hierarchical' requires institution labels")
        if len(groups) != n:
            raise InvalidOpinion(f"expected {n} group labels, got {len(groups)}")
        sizes: Dict[object, int] = {}
        for g in groups:
            sizes[g] = sizes.get(g, 0) + 1
        k = len(sizes)
        # Each institution carries total weight 1/K, split among its sources.
        src = [1.0 / (k * sizes[g]) for g in groups]
        raw = {(i, j): src[i] * src[j] for i, j in pairs}

    else:
        raise InvalidOpinion(f"unknown combine mode: {combine!r}")

    total = sum(raw.values())
    if total <= 0:
        raise InvalidOpinion("pair weights sum to zero; check weights or groups")
    return {p: w / total for p, w in raw.items()}


# =====================================================
# 4) AGGREGATION
# =====================================================

def aggregate(
    opinions: Sequence[IFV],
    shots: int = 20000,
    exact: bool = False,
    combine: str = "uniform",
    seed: Optional[int] = None,
    weights: Optional[Sequence[float]] = None,
    groups: Optional[Sequence[object]] = None,
) -> Tuple[IFV, List[Dict[str, object]], Dict[str, object]]:
    """
    Returns (C~, steps, stats).

    steps records every pair -- indices, weight and agreement degrees -- for
    explainability. stats carries the 95% intervals and, when groups are
    given, the intra/inter-institutional decomposition.

    Runs C(n,2) circuits, each with its own random stream derived from
    `seed`. There is no special path by n: the pipeline is the same for every
    n >= 2, and n = 2 is the case P = C(2,2) = 1 in which the aggregate
    reduces to the single boxplus_I(x~_1, x~_2).
    """
    n = len(opinions)
    if n < 2:
        raise InvalidOpinion(
            f"at least two opinions are required; got {n}"
        )
    if combine not in COMBINE_MODES:
        raise InvalidOpinion(f"unknown combine mode: {combine!r}")

    if combine == "tnorm":
        return _aggregate_tnorm(opinions, shots, exact, seed)

    w = pair_weights(n, combine=combine, weights=weights, groups=groups)
    pairs = sorted(w)
    seeds = derive_seeds(seed, len(pairs)) if seed is not None else [None] * len(pairs)

    steps: List[Dict[str, object]] = []
    for (i, j), sd in zip(pairs, seeds):
        if exact:
            mu, nu = xnor_analytic(opinions[i], opinions[j])
            p11, eff = mu * nu, 0          # eff = 0 signals zero variance
        else:
            r = run_xnor_joint(opinions[i], opinions[j], shots=shots, seed=sd)
            mu, nu, p11, eff = r["mu"], r["nu"], r["p11"], r["shots"]

        steps.append({
            "pair_i": i + 1,
            "pair_j": j + 1,
            "weight": w[(i, j)],
            "mu": mu,
            "nu": nu,
            "pi": max(0.0, 1.0 - mu - nu),
            "same_group": None if groups is None else groups[i] == groups[j],
            "_p11": p11,
            "_shots": eff,
        })

    mu_c = sum(s["weight"] * s["mu"] for s in steps)
    nu_c = sum(s["weight"] * s["nu"] for s in steps)
    stats = _stats(steps, mu_c, nu_c, groups)

    for s in steps:                          # internal fields never leave the API
        s.pop("_p11"), s.pop("_shots")

    return (mu_c, nu_c), steps, stats


def _aggregate_tnorm(opinions, shots, exact, seed):
    """T_I conjunction over the pairs: requires ALL of them to agree. Conservative."""
    pairs = list(itertools.combinations(range(len(opinions)), 2))
    seeds = derive_seeds(seed, len(pairs)) if seed is not None else [None] * len(pairs)

    steps, acc = [], None
    for (i, j), sd in zip(pairs, seeds):
        if exact:
            mu, nu = xnor_analytic(opinions[i], opinions[j])
        else:
            r = run_xnor_joint(opinions[i], opinions[j], shots=shots, seed=sd)
            mu, nu = r["mu"], r["nu"]
        steps.append({
            "pair_i": i + 1, "pair_j": j + 1, "weight": 1.0 / len(pairs),
            "mu": mu, "nu": nu, "pi": max(0.0, 1.0 - mu - nu), "same_group": None,
        })
        acc = (mu, nu) if acc is None else T_I_P(acc, (mu, nu))

    return acc, steps, _empty_stats(acc)


def _empty_stats(c: IFV) -> Dict[str, object]:
    return {
        "mu_ci95": 0.0, "nu_ci95": 0.0, "pi_ci95": 0.0,
        "intra": None, "inter": None, "n_pairs": 0,
    }


def _stats(steps, mu_c, nu_c, groups) -> Dict[str, object]:
    """95% intervals and intra/inter decomposition."""
    var_mu = var_nu = cov = 0.0
    for s in steps:
        shots = s["_shots"]
        if not shots:                        # exact mode: no sampling variance
            continue
        w2 = s["weight"] ** 2
        mu, nu, p11 = s["mu"], s["nu"], s["_p11"]
        var_mu += w2 * mu * (1 - mu) / shots
        var_nu += w2 * nu * (1 - nu) / shots
        cov += w2 * (p11 - mu * nu) / shots

    var_pi = max(0.0, var_mu + var_nu + 2 * cov)

    def _mean(sel):
        chosen = [s for s in steps if s["same_group"] is sel]
        if not chosen:
            return None
        tot = sum(s["weight"] for s in chosen)
        if tot <= 0:
            return None
        return {
            "mu": sum(s["weight"] * s["mu"] for s in chosen) / tot,
            "nu": sum(s["weight"] * s["nu"] for s in chosen) / tot,
            "n_pairs": len(chosen),
            "weight": tot,
        }

    return {
        "mu_ci95": _Z95 * math.sqrt(var_mu),
        "nu_ci95": _Z95 * math.sqrt(var_nu),
        "pi_ci95": _Z95 * math.sqrt(var_pi),
        "intra": None if groups is None else _mean(True),
        "inter": None if groups is None else _mean(False),
        "n_pairs": len(steps),
    }


# =====================================================
# 5) CLINICAL DECISION
# =====================================================

def decide(
    opinions: Sequence[IFV],
    shots: int = 20000,
    exact: bool = False,
    combine: str = "uniform",
    seed: Optional[int] = None,
    weights: Optional[Sequence[float]] = None,
    groups: Optional[Sequence[object]] = None,
) -> Dict[str, object]:
    """
    Maps n opinions onto a clinical action. The three scores sum to 1 by
    construction. Ties between treat and do-not-treat resolve to further
    examinations, which is the conservative choice.
    """
    for i, (mu, nu) in enumerate(opinions, start=1):
        validate_opinion(mu, nu, label=f"opinion[{i}]")

    consensus, steps, stats = aggregate(
        opinions, shots=shots, exact=exact, combine=combine,
        seed=seed, weights=weights, groups=groups,
    )
    mu_c, nu_c = consensus
    pi_c = max(0.0, 1.0 - mu_c - nu_c)

    # The mean evidence uses the same per-source weights as the aggregation.
    src_w = None
    if combine == "weighted":
        src_w = list(weights)
    elif combine == "hierarchical" and groups is not None:
        sizes: Dict[object, int] = {}
        for g in groups:
            sizes[g] = sizes.get(g, 0) + 1
        k = len(sizes)
        src_w = [1.0 / (k * sizes[g]) for g in groups]

    evidence = mean_ifv(opinions, src_w)
    s_e = score(evidence)

    scores = {
        DECISION_TREAT: mu_c * (1.0 + s_e) / 2.0,
        DECISION_DO_NOT_TREAT: mu_c * (1.0 - s_e) / 2.0,
        DECISION_REQUEST_EXAMS: 1.0 - mu_c,
    }

    best = max(scores.values())
    winners = [k for k, v in scores.items() if best - v <= _TOL]
    label = DECISION_REQUEST_EXAMS if len(winners) > 1 else winners[0]

    return {
        "decision": label,
        "scores": scores,
        "consensus": {"mu": mu_c, "nu": nu_c, "pi": pi_c},
        "uncertainty": {
            "mu_ci95": stats["mu_ci95"],
            "nu_ci95": stats["nu_ci95"],
            "pi_ci95": stats["pi_ci95"],
        },
        "breakdown": {"intra": stats["intra"], "inter": stats["inter"]},
        "evidence": {
            "mu": evidence[0], "nu": evidence[1],
            "score": s_e, "accuracy": accuracy(evidence),
        },
        "steps": steps,
        "n_opinions": len(opinions),
        "combine": combine,
        "provenance": (
            {"backend": "closed form (Proposition 2)", "shots": None, "seed": None}
            if exact else provenance(shots, seed)
        ),
    }
