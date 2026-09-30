"""
FuzzyQ-MedAgg -- API REST.

Exposes the Evidence Fusion Core, which aggregates n intuitionistic fuzzy
opinions through the intuitionistic fuzzy XNOR operator (QIF-XNOR Engine) and
returns a clinical action.
"""

from typing import Dict, List, Literal, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field, model_validator

import fusion
from fusion import InvalidOpinion

# =====================================================
# 1) MODELS
# =====================================================

InternalDecisionLabel = Literal["DO_NOT_TREAT", "REQUEST_EXAMS", "TREAT"]
ExternalDecisionLabel = Literal["TREAT", "DO NOT TREAT", "REQUEST EXAMS"]

_EXTERNAL = {
    "TREAT": "TREAT",
    "DO_NOT_TREAT": "DO NOT TREAT",
    "REQUEST_EXAMS": "REQUEST EXAMS",
}


class Opinion(BaseModel):
    """Intuitionistic fuzzy value (mu, nu) reported by a diagnostic source."""

    mu: float = Field(..., ge=0.0, le=1.0, description="membership degree")
    nu: float = Field(..., ge=0.0, le=1.0, description="non-membership degree")
    label: Optional[str] = None
    source_id: Optional[str] = Field(None, description="identifier of the source")
    institution_id: Optional[str] = Field(
        None, description="institution; used by combine='hierarchical'"
    )
    weight: float = Field(1.0, ge=0.0, description="weight; used by combine='weighted'")

    @model_validator(mode="after")
    def _check_intuitionistic(self):
        if self.mu + self.nu > 1.0 + 1e-9:
            raise ValueError(
                f"mu+nu={self.mu + self.nu:.4f} > 1 violates the intuitionistic condition"
            )
        return self


class DecisionRequest(BaseModel):
    """
    Accepts a list of n opinions or, for backward compatibility, the
    doctor1/doctor2 pair used by the previous version of the interface.
    """

    opinions: Optional[List[Opinion]] = None
    doctor1: Optional[Opinion] = None
    doctor2: Optional[Opinion] = None
    shots: int = Field(20000, ge=100, le=1_000_000)
    combine: Literal["uniform", "weighted", "hierarchical", "tnorm"] = "uniform"
    exact: bool = Field(
        False,
        description=(
            "Evaluate the closed form of Proposition 2 instead of sampling the "
            "circuit. Deterministic and free of sampling error; the provenance "
            "block reports which path was taken."
        ),
    )
    seed: Optional[int] = Field(
        None,
        ge=0,
        le=2**31 - 2,
        description=(
            "Omitting it keeps the run stochastic. Providing it makes the "
            "result reproducible bit for bit, with one derived stream per pair."
        ),
    )

    @model_validator(mode="after")
    def _normalize(self):
        if self.opinions is None:
            legacy = [o for o in (self.doctor1, self.doctor2) if o is not None]
            if not legacy:
                raise ValueError("provide either 'opinions' or 'doctor1'/'doctor2'")
            self.opinions = legacy
        if len(self.opinions) < 2:
            raise ValueError(
                f"at least two opinions are required; got {len(self.opinions)}"
            )
        return self


class IFVOut(BaseModel):
    mu: float
    nu: float
    pi: float


class ProvenanceOut(BaseModel):
    """Makes a published number traceable back to the run that produced it."""

    qiskit: Optional[str] = None
    qiskit_aer: Optional[str] = None
    backend: str
    shots: Optional[int] = None
    seed: Optional[int] = None


class EvidenceOut(BaseModel):
    mu: float
    nu: float
    score: float
    accuracy: float


class PairAgreementOut(BaseModel):
    """boxplus_I agreement between one pair of opinions (1-based indices)."""

    pair_i: int
    pair_j: int
    weight: float
    mu: float
    nu: float
    pi: float
    same_group: Optional[bool] = None


class UncertaintyOut(BaseModel):
    """
    Half-width of the 95% interval of each aggregated degree.

    mu_C and nu_C are weighted means of binomial proportions estimated from
    `shots` samples per pair; pi_C = 1 - mu_C - nu_C includes the covariance
    between mu and nu, which are read from the same measurement.
    """

    mu_ci95: float
    nu_ci95: float
    pi_ci95: float


class GroupAgreementOut(BaseModel):
    mu: float
    nu: float
    n_pairs: int
    weight: float


class BreakdownOut(BaseModel):
    """Agreement split into same-institution and cross-institution pairs."""

    intra: Optional[GroupAgreementOut] = None
    inter: Optional[GroupAgreementOut] = None


class DecisionResponse(BaseModel):
    decision: ExternalDecisionLabel
    scores: Dict[str, float]
    consensus: IFVOut
    uncertainty: UncertaintyOut
    breakdown: BreakdownOut
    evidence: EvidenceOut
    steps: List[PairAgreementOut]
    n_opinions: int
    combine: str
    provenance: ProvenanceOut
    interpretation: str

    # Flat fields kept for the existing interface.
    mu_xnor: float
    nu_xnor: float
    pi_xnor: float


# =====================================================
# 2) APP
# =====================================================

app = FastAPI(
    title="FuzzyQ-MedAgg",
    description=(
        "Agreement-based aggregation of intuitionistic fuzzy opinions "
        "through an XNOR operator realised as a quantum circuit."
    ),
    version="2.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["POST", "GET"],
    allow_headers=["*"],
)


def _interpretation(result: dict) -> str:
    c = result["consensus"]
    e = result["evidence"]
    return (
        f"Agreement mu={c['mu']:.2f}, conflict nu={c['nu']:.2f}, "
        f"hesitation pi={c['pi']:.2f} (95% CI +/-{result['uncertainty']['mu_ci95']:.3f} on mu) "
        f"across {result['n_opinions']} sources "
        f"({len(result['steps'])} pairwise comparisons, {result['combine']} weighting). "
        f"Evidence score s={e['score']:+.2f} "
        f"({'supports' if e['score'] > 0 else 'opposes' if e['score'] < 0 else 'neutral on'} "
        f"the diagnosis). Final decision: {_EXTERNAL[result['decision']]}."
    )


@app.get("/health")
def health() -> dict:
    import qif_xnor

    return {
        "status": "ok",
        "version": app.version,
        **qif_xnor.provenance(shots=0),
    }


@app.post("/decide", response_model=DecisionResponse)
def decide(payload: DecisionRequest) -> DecisionResponse:
    opinions = [(o.mu, o.nu) for o in payload.opinions]
    weights = [o.weight for o in payload.opinions]
    groups = [
        o.institution_id if o.institution_id is not None else f"__source_{i}"
        for i, o in enumerate(payload.opinions)
    ]

    try:
        result = fusion.decide(
            opinions,
            shots=payload.shots,
            combine=payload.combine,
            seed=payload.seed,
            exact=payload.exact,
            weights=weights,
            groups=groups,
        )
    except InvalidOpinion as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    c = result["consensus"]
    return DecisionResponse(
        decision=_EXTERNAL[result["decision"]],
        scores={
            "treat": result["scores"]["TREAT"],
            "doNotTreat": result["scores"]["DO_NOT_TREAT"],
            "requestExams": result["scores"]["REQUEST_EXAMS"],
        },
        consensus=IFVOut(**c),
        uncertainty=UncertaintyOut(**result["uncertainty"]),
        breakdown=BreakdownOut(**result["breakdown"]),
        evidence=EvidenceOut(**result["evidence"]),
        steps=[PairAgreementOut(**s) for s in result["steps"]],
        n_opinions=result["n_opinions"],
        combine=result["combine"],
        provenance=ProvenanceOut(**result["provenance"]),
        interpretation=_interpretation(result),
        mu_xnor=c["mu"],
        nu_xnor=c["nu"],
        pi_xnor=c["pi"],
    )
