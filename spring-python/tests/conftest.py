"""Shared fixtures for the SPRING port test-suite."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

GOLDEN_DIR = Path(os.environ.get("SPRING_GOLDEN_DIR", Path(__file__).parent / "golden"))


def _case_params(name: str) -> dict:
    return json.loads((GOLDEN_DIR / name / "params.json").read_text())


_ALL_CASES = sorted(
    path.name
    for path in GOLDEN_DIR.iterdir()
    if path.is_dir() and (path / "params.json").is_file()
)
# Full-pipeline cases (SPRING + StARS). The dedicated `kind="latentcor"` fixtures
# compare only the latent-correlation stage, because `Rmethod="original"` is far
# too slow to run a full StARS golden on.
CASES = [name for name in _ALL_CASES if _case_params(name).get("kind") != "latentcor"]
LATENTCOR_CASES = [
    name for name in _ALL_CASES if _case_params(name).get("kind") == "latentcor"
]


def _read_matrix(path: Path) -> np.ndarray:
    """Read a whitespace/comma-separated full-precision matrix."""
    data = np.loadtxt(path, delimiter=",", ndmin=2)
    return data


def load_case(name: str) -> dict:
    """Load all golden artefacts of one case."""
    d = GOLDEN_DIR / name
    if not d.is_dir():
        pytest.skip(f"golden case {name!r} not available")
    params = json.loads((d / "params.json").read_text())
    out = {
        "dir": d,
        "params": params,
        "input": _read_matrix(d / "input.csv"),
        "qdat": _read_matrix(d / "qdat.csv"),
        "lambdaseq": np.loadtxt(d / "lambdaseq.csv", delimiter=",").ravel(),
        "K": _read_matrix(d / "full_K.csv"),
        "Rpointwise": _read_matrix(d / "full_Rpointwise.csv"),
        "R": _read_matrix(d / "full_R.csv"),
        "summary": np.loadtxt(d / "summary.csv", delimiter=",").ravel(),
        "opt_index": int((d / "opt_index.txt").read_text().strip()),
        "subsample_indices": _read_matrix(d / "subsample_indices.csv").astype(int) - 1,
        "refit": _read_matrix(d / "adjacency_refit.csv"),
    }
    nlam = len(out["lambdaseq"])
    out["merge"] = [_read_matrix(d / f"merge_{i:02d}.csv") for i in range(1, nlam + 1)]
    out["beta_full"] = [_read_matrix(d / f"beta_full_{i:02d}.csv") for i in range(1, nlam + 1)]
    out["adjacency_full"] = [
        _read_matrix(d / f"adjacency_full_{i:02d}.csv") for i in range(1, nlam + 1)
    ]
    out["subsample_R"] = [
        _read_matrix(d / f"subsample_R_{k:02d}.csv")
        for k in range(1, len(out["subsample_indices"]) + 1)
    ]
    return out


def load_latentcor_case(name: str) -> dict:
    """Load a dedicated latent-correlation-only golden case."""
    d = GOLDEN_DIR / name
    if not d.is_dir():
        pytest.skip(f"latentcor golden case {name!r} not available")
    return {
        "dir": d,
        "params": json.loads((d / "params.json").read_text()),
        "input": _read_matrix(d / "input.csv"),
        "qdat": _read_matrix(d / "qdat.csv"),
        "K": _read_matrix(d / "full_K.csv"),
        "Rpointwise": _read_matrix(d / "full_Rpointwise.csv"),
        "R": _read_matrix(d / "full_R.csv"),
    }


@pytest.fixture(params=CASES)
def case(request):
    """Parametrised fixture over the golden cases."""
    return load_case(request.param)
