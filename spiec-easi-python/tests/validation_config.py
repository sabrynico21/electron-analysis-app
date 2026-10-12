from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DeterministicTolerance:
    """Numerical tolerance for deterministic component-level checks."""

    rtol: float
    atol: float


@dataclass(frozen=True)
class StatisticalBaselineTolerance:
    """Acceptance margins for stochastic baseline comparisons."""

    edge_jaccard_distance_margin: float
    degree_mae_margin: float
    lambda_rank_margin: float


# Fase 1 (component-level lasso parity): coefficients from sklearn-based
# `_fit_single_neighborhood` vs R glmnet with explicit lambda, strict optimizer
# tolerance and fixed solver settings. These are declared a priori and are not
# tuned post-hoc to force passing tests.
LASSO_COEFFICIENT_TOLERANCE = DeterministicTolerance(
    rtol=1e-4,
    atol=1e-6,
)


# Coefficients smaller than this are treated as structural zeros when evaluating
# support exactness; chosen to be close to numerical precision while avoiding
# unstable sign flips around machine epsilon.
SUPPORT_ZERO_THRESHOLD = 1e-10


# Relative errors are unstable for near-zero coefficients; this floor defines
# "informative" non-zero coefficients for rtol checks, while atol remains
# enforced on all coefficients.
RELATIVE_ERROR_MIN_MAGNITUDE = 1e-3


# Fase 1 coverage requirements.
MIN_LASSO_TARGET_COMPARISONS = 20
MIN_LASSO_LAMBDAS_PER_TARGET = 5


# Fase 2 (fixed-lambda MB network parity): any adjacency discordance is only
# acceptable for edges whose symmetrized weights are numerically negligible.
FIXED_LAMBDA_AMBIGUOUS_EDGE_WEIGHT_MAX = 1e-6
FIXED_LAMBDA_WEIGHT_TOLERANCE = DeterministicTolerance(
    rtol=1e-4,
    atol=1e-6,
)


# Fase 3 (shared subsamples): with identical CLR matrix, lambda path, and
# subsamples, StARS instability summaries should match nearly exactly across
# Python and R implementations.
STARS_SHARED_SUBSAMPLES_TOLERANCE = DeterministicTolerance(
    rtol=1e-6,
    atol=1e-8,
)


# Fase 3 scenario coverage.
MIN_SHARED_SUBSAMPLE_SCENARIOS = 5


# Fase 4 (stochastic baseline): cross-stack variability (Py-vs-R with
# independent subsampling RNG) is compared against within-stack variability
# (Py-vs-Py and R-vs-R). Margins are declared a priori and not tuned on outcomes.
STOCHASTIC_BASELINE_TOLERANCE = StatisticalBaselineTolerance(
    edge_jaccard_distance_margin=0.08,
    degree_mae_margin=0.20,
    lambda_rank_margin=0.10,
)


# Fase 4 minimum benchmark coverage.
MIN_STOCHASTIC_BASELINE_SCENARIOS = 3
MIN_STOCHASTIC_BASELINE_SEEDS = 4


# Fase 5 (expanded synthetic benchmark) minimum topology coverage.
MIN_EXPANDED_SCENARIO_FAMILIES = 3
MIN_EXPANDED_SCENARIOS = 6
