/**
 * Client for the FuzzyQ-MedAgg backend.
 *
 * The base URL comes from NEXT_PUBLIC_API_URL so the app can be deployed
 * without editing source. See .env.example.
 */

export const API_URL =
  process.env.NEXT_PUBLIC_API_URL ?? "http://127.0.0.1:8000"

/** Intuitionistic fuzzy value produced by one diagnostic source. */
export interface Opinion {
  id: string
  label: string
  mu: number
  nu: number
  /** Site the source belongs to. Drives the "by institution" weighting. */
  institution?: string
}

/**
 * uniform      every pair counts equally
 * hierarchical every institution carries the same total weight, whatever its size
 * weighted     per-source weights (API only; not exposed in this interface)
 */
export type CombineMode = "uniform" | "hierarchical" | "weighted"

export type Decision = "TREAT" | "DO NOT TREAT" | "REQUEST EXAMS"

/** Agreement between one pair of opinions (1-based indices). */
export interface PairAgreement {
  pair_i: number
  pair_j: number
  weight: number
  mu: number
  nu: number
  pi: number
  same_group: boolean | null
}

/** Makes a published number traceable back to the run that produced it. */
export interface Provenance {
  qiskit?: string | null
  qiskit_aer?: string | null
  backend: string
  shots?: number | null
  seed?: number | null
}

/** Half-width of the 95% interval on each aggregated degree. */
export interface Uncertainty {
  mu_ci95: number
  nu_ci95: number
  pi_ci95: number
}

export interface GroupAgreement {
  mu: number
  nu: number
  n_pairs: number
  weight: number
}

/** Agreement split into same-institution and cross-institution pairs. */
export interface Breakdown {
  intra: GroupAgreement | null
  inter: GroupAgreement | null
}

export interface DecisionResult {
  decision: Decision
  scores: { treat: number; doNotTreat: number; requestExams: number }
  consensus: { mu: number; nu: number; pi: number }
  uncertainty: Uncertainty
  breakdown: Breakdown
  evidence: { mu: number; nu: number; score: number; accuracy: number }
  steps: PairAgreement[]
  n_opinions: number
  combine: CombineMode
  provenance: Provenance
  interpretation: string
}

/** An opinion is valid only inside U~ = {(mu, nu) : mu + nu <= 1}. */
export function isValid(o: Pick<Opinion, "mu" | "nu">): boolean {
  return o.mu >= 0 && o.nu >= 0 && o.mu <= 1 && o.nu <= 1 && o.mu + o.nu <= 1
}

export function hesitation(o: Pick<Opinion, "mu" | "nu">): number {
  return Math.max(0, 1 - o.mu - o.nu)
}

/** Score function of the Xu-Yager order: s(x) = mu - nu, in [-1, 1]. */
export function score(o: Pick<Opinion, "mu" | "nu">): number {
  return o.mu - o.nu
}

/** Accuracy function of the Xu-Yager order: h(x) = mu + nu, in [0, 1]. */
export function accuracy(o: Pick<Opinion, "mu" | "nu">): number {
  return o.mu + o.nu
}

/** Distinct institutions named across the given sources. */
export function institutionsOf(opinions: Opinion[]): string[] {
  const seen = new Set<string>()
  for (const o of opinions) {
    const name = o.institution?.trim()
    if (name) seen.add(name)
  }
  return [...seen].sort()
}

export async function requestDecision(
  opinions: Opinion[],
  shots = 20000,
  seed?: number,
  combine: CombineMode = "uniform",
  /** Closed form of Proposition 2 instead of sampling; deterministic. */
  exact = false,
): Promise<DecisionResult> {
  const res = await fetch(`${API_URL}/decide`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      opinions: opinions.map(({ mu, nu, label, institution }) => ({
        mu,
        nu,
        label,
        institution_id: institution?.trim() ? institution.trim() : null,
      })),
      shots,
      combine,
      exact,
      ...(seed === undefined ? {} : { seed }),
    }),
  })

  if (!res.ok) {
    let detail = `Request failed with status ${res.status}`
    try {
      const body = await res.json()
      if (typeof body.detail === "string") detail = body.detail
      else if (Array.isArray(body.detail) && body.detail[0]?.msg)
        detail = body.detail[0].msg
    } catch {
      // response had no JSON body; keep the status-based message
    }
    throw new Error(detail)
  }

  return res.json()
}
