"""
Shared infrastructure for the experimental evaluation.

This package does NOT alter production behaviour: it imports qif_xnor, fusion
and resources and only observes them.
"""

import csv
import os
import platform
import sys
from typing import Dict, List, Sequence, Tuple

_HERE = os.path.dirname(os.path.abspath(__file__))
_BACKEND = os.path.dirname(os.path.dirname(_HERE))
if _BACKEND not in sys.path:
    sys.path.insert(0, _BACKEND)

OUT_DIR = _HERE

# ---------------------------------------------------------------
# Protocol: shot budgets and seeds
# ---------------------------------------------------------------
SHOT_BUDGETS = (2_000, 8_000, 32_000)
PRODUCTION_DEFAULT_SHOTS = 20_000

SEEDS = {
    "analytical_vs_circuit": 20260901,   # base seed per case, offset by the index
    "coverage": 70000,                   # base seed per repetition
    "multi_source": 4242,
    "runtime": 909,
    "opinion_sampling": 13337,           # generation of random valid opinions
}

# Grid resolutions over U~
GRID_STEP_PROPERTIES = 0.05   # 1/20 -> 231 pontos validos
GRID_STEP_CLOSURE = 0.02      # 1/50 -> 1326 pontos validos

TOL_EXACT = 1e-12             # floating-point equality comparisons
TOL_SUM = 1e-9                # slack on the condition mu + nu <= 1


# ---------------------------------------------------------------
# Representative cases (actual tuples, no invented labels)
# ---------------------------------------------------------------
# regime, x~ = (mu_x, nu_x), y~ = (mu_y, nu_y)
REPRESENTATIVE_CASES: List[Tuple[str, Tuple[float, float], Tuple[float, float]]] = [
    ("support-oriented agreement",   (0.90, 0.05), (0.85, 0.10)),
    ("rejection-oriented agreement", (0.05, 0.90), (0.10, 0.85)),
    ("strong disagreement",          (0.90, 0.05), (0.05, 0.90)),
    ("high hesitation",              (0.30, 0.30), (0.30, 0.30)),
    ("asymmetric evidence",          (0.80, 0.15), (0.35, 0.45)),
    ("zero-hesitation boundary",     (0.60, 0.40), (0.55, 0.45)),
    ("near-vertex",                  (0.99, 0.01), (0.98, 0.02)),
]


def valid_grid(step: float) -> List[Tuple[float, float]]:
    """Points (mu, nu) of U~ on a regular grid of step `step`."""
    k = int(round(1.0 / step))
    vals = [i / k for i in range(k + 1)]
    return [(a, b) for a in vals for b in vals if a + b <= 1.0 + TOL_SUM]


def sample_opinions(n: int, seed: int) -> List[Tuple[float, float]]:
    """
    n valid opinions sampled uniformly from the simplex of U~.

    Uses a dedicated Generator, without touching numpy's global state.
    """
    import numpy as np

    rng = np.random.default_rng(seed)
    out = []
    while len(out) < n:
        mu, nu = rng.random(2)
        if mu + nu <= 1.0:
            out.append((float(mu), float(nu)))
    return out


def write_csv(name: str, fieldnames: Sequence[str], rows: List[Dict]) -> str:
    path = os.path.join(OUT_DIR, name)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(fieldnames))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


def environment() -> Dict[str, object]:
    """Environment relevant to reproducibility, read from the system."""
    import fastapi
    import numpy
    import pydantic
    import qiskit
    import qiskit_aer

    from qif_xnor import BACKEND_NAME

    return {
        "python": sys.version.split()[0],
        "python_implementation": platform.python_implementation(),
        "qiskit": qiskit.__version__,
        "qiskit_aer": qiskit_aer.__version__,
        "numpy": numpy.__version__,
        "fastapi": fastapi.__version__,
        "pydantic": str(pydantic.VERSION),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "backend": BACKEND_NAME,
    }


def banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(f"  {title}")
    print("=" * 78)
