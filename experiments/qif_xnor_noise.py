"""QIF-XNOR depolarizing-noise study, v2.
(a) Exact noisy marginals via density-matrix simulation (isolates gate-noise effect, no sampling).
(b) 30 independent shot-based repetitions per noise level (20,000 shots, simulator seeds 1..30):
    MAE reported as mean +- std over repetitions.
Circuit built from the paper's description; transpiled to {cx,u}, optimization_level=1, seed_transpiler=1.
Writes noise_results.csv and prints a summary.  Versions are printed for the record."""
import numpy as np, csv, itertools, platform
import qiskit, qiskit_aer
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error

def SP(a, b): return a + b - a * b
def xnor(x, y):
    (x1, x2), (y1, y2) = x, y
    return SP(x2, y1) * SP(x1, y2), SP(x1 * y2, x2 * y1)
def sblock(qc, i, j, k):
    qc.x(i); qc.x(j); qc.ccx(i, j, k); qc.x(k); qc.x(i); qc.x(j)
def build(x, y, measure=True):
    qc = QuantumCircuit(10, 2 if measure else 0)
    for q, p in zip(range(4), [x[0], x[1], y[0], y[1]]):
        qc.ry(2 * np.arcsin(np.sqrt(p)), q)
    sblock(qc, 1, 2, 4); sblock(qc, 0, 3, 5); qc.ccx(4, 5, 8)
    qc.ccx(0, 3, 6); qc.ccx(1, 2, 7); sblock(qc, 6, 7, 9)
    if measure: qc.measure(8, 0); qc.measure(9, 1)
    return qc

PAIRS = [((0.85, 0.10), (0.80, 0.02)), ((0.10, 0.85), (0.02, 0.80)), ((0.85, 0.10), (0.10, 0.85)),
         ((0.30, 0.20), (0.25, 0.25)), ((0.90, 0.05), (0.60, 0.20)), ((1.0, 0.0), (0.0, 1.0)), ((0.5, 0.5), (0.5, 0.5))]
LEVELS = [0.0, 1e-4, 5e-4, 1e-3, 5e-3]
SHOTS, REPS = 20000, 30

def noise_model(p1):
    if p1 == 0: return None
    nm = NoiseModel()
    nm.add_all_qubit_quantum_error(depolarizing_error(p1, 1), ["u"])
    nm.add_all_qubit_quantum_error(depolarizing_error(10 * p1, 2), ["cx"])
    return nm

def marginals_from_probs(probs):
    # probs over 10-qubit basis states, little-endian index: bit k of index = qubit k
    idx = np.arange(len(probs))
    mu = probs[(idx >> 8) & 1 == 1].sum(); nu = probs[(idx >> 9) & 1 == 1].sum()
    return mu, nu

def exact(p1):
    sim = AerSimulator(method="density_matrix", noise_model=noise_model(p1))
    errs = []
    for x, y in PAIRS:
        qc = transpile(build(x, y, measure=False), basis_gates=["cx", "u"], optimization_level=1, seed_transpiler=1)
        qc.save_probabilities()
        probs = np.asarray(sim.run(qc).result().data(0)["probabilities"])
        mu, nu = marginals_from_probs(probs); m, n = xnor(x, y)
        errs.append((abs(mu - m), abs(nu - n)))
    e = np.array(errs); return e[:, 0].mean(), e[:, 1].mean()

def sampled(p1, seed):
    sim = AerSimulator(noise_model=noise_model(p1), seed_simulator=seed)
    errs = []
    for x, y in PAIRS:
        qc = transpile(build(x, y), basis_gates=["cx", "u"], optimization_level=1, seed_transpiler=1)
        c = sim.run(qc, shots=SHOTS).result().get_counts()
        mu = sum(v for k, v in c.items() if k[-1] == "1") / SHOTS
        nu = sum(v for k, v in c.items() if k[-2] == "1") / SHOTS
        m, n = xnor(x, y); errs.append((abs(mu - m), abs(nu - n)))
    e = np.array(errs); return e[:, 0].mean(), e[:, 1].mean()

print(f"python {platform.python_version()}  qiskit {qiskit.__version__}  qiskit-aer {qiskit_aer.__version__}")
rows = []
for p1 in LEVELS:
    ex_mu, ex_nu = exact(p1)
    reps = np.array([sampled(p1, s) for s in range(1, REPS + 1)])
    row = dict(p1=p1, p2=10 * p1, exact_mae_mu=ex_mu, exact_mae_nu=ex_nu,
               samp_mae_mu_mean=reps[:, 0].mean(), samp_mae_mu_std=reps[:, 0].std(ddof=1),
               samp_mae_nu_mean=reps[:, 1].mean(), samp_mae_nu_std=reps[:, 1].std(ddof=1))
    rows.append(row)
    print(f"p1={p1:.0e}: exact MAE mu={ex_mu:.4f} nu={ex_nu:.4f} | sampled ({REPS}x{SHOTS}) "
          f"mu={row['samp_mae_mu_mean']:.4f}±{row['samp_mae_mu_std']:.4f} nu={row['samp_mae_nu_mean']:.4f}±{row['samp_mae_nu_std']:.4f}")
with open("noise_results.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
