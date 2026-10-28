"""StARS stability selection (port of ``pulsar`` 0.3.13 for SPRING's usage).

Only the ``criterion = "stars"`` path (no B-StARS bounds) is implemented, which
is exactly what SPRING uses. Randomness is confined to the subsample indices,
which can be injected (``subsample_indices``) so results become reproducible
across R and Python.

R sources:
    pulsar/R/mcPulsarSelect.R  (pulsar)
    pulsar/R/utilities.R       (.ratcheck, .starsind, .tlist)
    pulsar/R/refit.R           (refit.pulsar)
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "ratcheck",
    "starsind",
    "stars_stability",
    "pulsar",
    "PulsarResult",
    "RefitResult",
]


def ratcheck(subsample_ratio: float | None, n: int) -> float:
    """``pulsar:::.ratcheck``: 0.8 if ``n <= 144`` else ``10*sqrt(n)/n``."""
    if subsample_ratio is None:
        if n > 144:
            return 10.0 * np.sqrt(n) / n
        return 0.8
    return float(subsample_ratio)


def starsind(summary, thresh: float, offset: int = 1) -> int:
    """``pulsar:::.starsind`` (returns a 1-based index)."""
    summary = np.asarray(summary, dtype=float)
    idx = np.flatnonzero(summary >= thresh)
    if idx.size:
        return max(int(idx[0]) + 1 - offset, 1)
    warnings.warn("Optimal lambda may be outside the supplied values", stacklevel=2)
    return int(summary.size)


@dataclass
class StarsResult:
    """StARS selection output (mirrors ``pulsar`` ``$stars``)."""

    merge: list[np.ndarray]
    summary: np.ndarray
    opt_index: int
    criterion: str = "stars.stability"
    thresh: float = 0.1


def stars_stability(premerge_reord, thresh: float, rep_num: int, p: int) -> StarsResult:
    """``pulsar:::stars.stability``.

    Parameters
    ----------
    premerge_reord : list of list
        ``[lambda][replicate]`` adjacency matrices.
    """
    merge = []
    for lam_paths in premerge_reord:
        acc = np.zeros_like(np.asarray(lam_paths[0], dtype=float))
        for path in lam_paths:
            acc += np.asarray(path, dtype=float)
        merge.append(acc / rep_num)

    summary = np.zeros(len(merge))
    for i, m in enumerate(merge):
        summary[i] = 4.0 * np.sum(m * (1.0 - m)) / (p * (p - 1))
    summary = np.maximum.accumulate(summary)
    opt_index = starsind(summary, thresh)
    return StarsResult(merge=merge, summary=summary, opt_index=opt_index, thresh=thresh)


@dataclass
class RefitResult:
    """Mirrors ``pulsar.refit`` output."""

    est: object
    refit: dict = field(default_factory=dict)


@dataclass
class PulsarResult:
    """Mirrors the ``pulsar`` object returned by SPRING."""

    stars: StarsResult
    est: object
    call: dict = field(default_factory=dict)


def pulsar(
    data: np.ndarray,
    fun: Callable,
    fargs: dict,
    criterion: str = "stars",
    thresh: float = 0.1,
    subsample_ratio: float | None = None,
    rep_num: int = 20,
    seed: int | None = None,
    ncores: int = 1,
    refit: bool = True,
    *,
    subsample_indices=None,
    rng: np.random.Generator | None = None,
) -> PulsarResult:
    """Run StARS stability selection over a regularization path.

    ``ncores`` is accepted for API compatibility; execution is serial (the
    estimator is deterministic, so parallelism cannot change the result).
    """
    data = np.asarray(data, dtype=float)
    n, p = data.shape
    if "lambda" not in fargs:
        raise ValueError("Error: missing members in fargs: lambda")
    lam = np.asarray(fargs["lambda"], dtype=float)
    if not np.all(lam == np.minimum.accumulate(lam)):
        warnings.warn(
            "Are you sure you don't want the lambda path to be monotonically decreasing",
            stacklevel=2,
        )
    if lam.size < 2:
        warnings.warn(
            "Only 1 value of lambda is given. Are you sure you want to do model selection?",
            stacklevel=2,
        )
    nlams = lam.size
    if criterion != "stars":
        raise NotImplementedError("Only criterion='stars' is supported.")

    ratio = ratcheck(subsample_ratio, n)
    size = int(np.floor(n * ratio))

    if subsample_indices is not None:
        ind_sample = [np.asarray(ix, dtype=int) for ix in subsample_indices]
        if len(ind_sample) != rep_num:
            raise ValueError(
                f"subsample_indices must contain rep.num={rep_num} index arrays."
            )
    else:
        if rng is None:
            rng = np.random.default_rng(seed)
        ind_sample = [rng.choice(n, size=size, replace=False) for _ in range(rep_num)]

    full_index = np.arange(n)
    seq = ([full_index] + ind_sample) if refit else ind_sample

    premerge = []
    full_model = None
    for k, idx in enumerate(seq):
        result = fun(data[idx], **fargs)
        if refit and k == 0:
            full_model = result
            continue
        premerge.append(result["path"])

    # reorder to [lambda][replicate]
    premerge_reord = [
        [premerge[j][i] for j in range(len(premerge))] for i in range(nlams)
    ]

    stars = stars_stability(premerge_reord, thresh, rep_num, p)
    if stars.opt_index == 1:
        direction = "larger" if np.any(stars.summary >= 0.1) else "smaller"
        warnings.warn(
            f"Optimal lambda may be {direction} than the supplied values", stacklevel=2
        )
    return PulsarResult(stars=stars, est=full_model)


def refit(obj: PulsarResult, criterion: str = "stars") -> RefitResult:
    """``pulsar::refit`` for the ``stars`` criterion."""
    est = obj.est
    out = RefitResult(est=est)
    optind = obj.stars.opt_index
    out.refit[criterion] = est["path"][optind - 1] if optind is not None else None
    return out
