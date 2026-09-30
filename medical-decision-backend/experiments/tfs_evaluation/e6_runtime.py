"""
E6 -- Computational resources and scalability.

Reports the cost of ONE pairwise circuit (fixed width and depth, since it
always evaluates exactly two opinions) and the growth of the number of
circuits with n. Times are measured separating classical fusion from quantum
simulation.
"""

import statistics
import time
from collections import Counter
from typing import Dict, List

from _common import PRODUCTION_DEFAULT_SHOTS, SEEDS, banner, write_csv

from qiskit import QuantumCircuit, transpile
from qiskit_aer import Aer

import resources as R
from fusion import aggregate, pair_weights
from qif_xnor import BACKEND_NAME, build_xnor_circuit, figure_ops

POOL = [(0.90, 0.05), (0.85, 0.10), (0.20, 0.70), (0.55, 0.35),
        (0.40, 0.50), (0.70, 0.15), (0.30, 0.60), (0.62, 0.28),
        (0.15, 0.75), (0.48, 0.42), (0.88, 0.08), (0.25, 0.65)]

N_VALUES = (2, 3, 4, 5, 8, 10, 12)
TIMING_REPS = 5


def main():
    banner("E6  COMPUTATIONAL RESOURCES AND SCALABILITY")

    qc = build_xnor_circuit((0.70, 0.20), (0.65, 0.25))
    fig = figure_ops(qc)
    full = Counter(ci.operation.name for ci in qc.data)
    logical = Counter(ci.operation.name for ci in fig)
    only = QuantumCircuit(qc.num_qubits, qc.num_clbits)
    for ci in fig:
        only.append(ci.operation, ci.qubits, ci.clbits)
    tc = transpile(qc, Aer.get_backend(BACKEND_NAME))

    print("  ONE pairwise circuit (width and depth independent of n):")
    print(f"    qubits                {qc.num_qubits}")
    print(f"    state preparation     {full['u']} u")
    print(f"    Pauli-X               {logical['x']}")
    print(f"    Toffoli (CCX)         {logical['ccx']}")
    print(f"    measurements          {full['measure']}")
    print(f"    total operations      {len(qc.data)}")
    print(f"    logical depth         {only.depth()}  (block T1--T11)")
    print(f"    full depth            {qc.depth()}")
    print(f"    transpilada (backend) {tc.depth()}  {dict(tc.count_ops())}")

    print("\n  estimates under explicit topologies (resources.py):")
    print(f"    {'topology':12}{'device qubits':>15}{'depth':>8}{'CNOTs':>8}")
    topo_rows = []
    for e in R.estimate_all(qc):
        dev = e["device_qubits"] if e["device_qubits"] else "free"
        print(f"    {e['topology']:12}{str(dev):>15}{e['depth']:>8}{e['cx']:>8}")
        topo_rows.append(e)

    rows: List[Dict] = []
    print(f"\n  scalability ({TIMING_REPS} repetitions per n,"
          f" {PRODUCTION_DEFAULT_SHOTS} shots):")
    print(f"    {'n':>3}{'P':>5}{'weights (ms)':>14}{'exact fusion':>15}"
          f"{'sampled total':>16}{'simulation':>13}{'ms/pair':>10}")
    for n in N_VALUES:
        ops = POOL[:n]
        P = n * (n - 1) // 2

        tw = []
        for _ in range(TIMING_REPS * 20):
            t = time.perf_counter()
            pair_weights(n, "uniform")
            tw.append(time.perf_counter() - t)
        t_w = statistics.median(tw)

        te = []
        for _ in range(TIMING_REPS * 5):
            t = time.perf_counter()
            aggregate(ops, exact=True)
            te.append(time.perf_counter() - t)
        t_exact = statistics.median(te)

        ts = []
        for k in range(TIMING_REPS):
            t = time.perf_counter()
            aggregate(ops, shots=PRODUCTION_DEFAULT_SHOTS,
                      seed=SEEDS["runtime"] + k)
            ts.append(time.perf_counter() - t)
        t_tot = statistics.median(ts)
        t_sim = max(0.0, t_tot - t_exact)

        rows.append({
            "n": n, "n_pairs": P, "shots": PRODUCTION_DEFAULT_SHOTS,
            "timing_repetitions": TIMING_REPS,
            "pair_weight_ms_median": f"{t_w * 1e3:.4f}",
            "classical_fusion_ms_median": f"{t_exact * 1e3:.3f}",
            "total_sampled_ms_median": f"{t_tot * 1e3:.1f}",
            "total_sampled_ms_mean": f"{statistics.fmean(ts) * 1e3:.1f}",
            "total_sampled_ms_stdev": (
                f"{statistics.stdev(ts) * 1e3:.1f}" if len(ts) > 1 else ""),
            "quantum_simulation_ms_median": f"{t_sim * 1e3:.1f}",
            "ms_per_pair": f"{t_sim * 1e3 / P:.1f}",
        })
        print(f"    {n:>3}{P:>5}{t_w*1e3:>13.4f}{t_exact*1e3:>13.3f}ms"
              f"{t_tot*1e3:>16.1f}ms{t_sim*1e3:>12.1f}ms{t_sim*1e3/P:>10.1f}")

    for r in rows:
        for e in topo_rows:
            r[f"circuit_depth_{e['topology'].replace('-', '_')}"] = e["depth"]
            r[f"circuit_cx_{e['topology'].replace('-', '_')}"] = e["cx"]
        r["circuit_qubits"] = qc.num_qubits
        r["circuit_x"] = logical["x"]
        r["circuit_ccx"] = logical["ccx"]
        r["circuit_prep_u"] = full["u"]
        r["circuit_measure"] = full["measure"]
        r["circuit_depth_logical_block"] = only.depth()
        r["circuit_depth_full"] = qc.depth()

    write_csv("runtime_scalability.csv", list(rows[0]), rows)
    print("\n  -> runtime_scalability.csv")
    print("  The number of circuits grows as P = n(n-1)/2 = O(n^2);"
          " each circuit has fixed width and depth.")
    return rows


if __name__ == "__main__":
    main()
