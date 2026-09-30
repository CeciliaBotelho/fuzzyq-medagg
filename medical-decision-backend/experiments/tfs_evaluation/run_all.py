"""
Runs the complete experimental evaluation in the order E1..E6.

    python3 run_all.py

No script alters production behaviour: all of them import qif_xnor, fusion and
resources and only observe them.
"""

import sys
import time

import e1_setup
import e2_properties
import e3_analytical_vs_circuit
import e4_coverage
import e5_multisource
import e6_runtime

STEPS = [
    ("E1 setup", e1_setup.main),
    ("E2 XNOR properties", e2_properties.main),
    ("E3 analytical vs circuit", e3_analytical_vs_circuit.main),
    ("E4 finite-shot coverage", e4_coverage.main),
    ("E5 multi-source aggregation", e5_multisource.main),
    ("E6 resources and scalability", e6_runtime.main),
]

if __name__ == "__main__":
    t0 = time.perf_counter()
    for name, fn in STEPS:
        t = time.perf_counter()
        fn()
        print(f"\n[{name} finished in {time.perf_counter() - t:.1f}s]",
              file=sys.stderr)
    print(f"\n[total {time.perf_counter() - t0:.1f}s]", file=sys.stderr)
