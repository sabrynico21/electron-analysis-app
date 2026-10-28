"""Command-line interface: JSON on stdin/``--input`` -> JSON on stdout.

Designed for the Electron/Node host (``child_process.spawn``). Nothing but JSON
is ever written to stdout; diagnostics go to stderr and a non-zero exit code is
returned on error.
"""

from __future__ import annotations

import argparse
import json
import sys

import numpy as np

from . import __version__
from .spring import spring
from .validation import validate_counts

__all__ = ["main"]

_ALIASES = {
    "lambda_min_ratio": "lambda_min_ratio",
    "lambda.min.ratio": "lambda_min_ratio",
    "subsample_ratio": "subsample_ratio",
    "subsample.ratio": "subsample_ratio",
    "rep_num": "rep_num",
    "rep.num": "rep_num",
    "Rmethod": "Rmethod",
    "use_nearPD": "use_nearPD",
    "Rtol": "Rtol",
}

_PARAMS = {
    "quantitative": False,
    "lambda_min_ratio": 1e-2,
    "nlambda": 20,
    "lambdaseq": None,
    "seed": 10010,
    "ncores": 1,
    "thresh": 0.1,
    "subsample_ratio": 0.8,
    "rep_num": 20,
    "Rtol": 1e-6,
    "verbose": False,
    "Rmethod": "approx",
    "use_nearPD": True,
    "nu": 0.001,
    "ratio": 0.9,
}


def _normalise_params(raw: dict) -> dict:
    params = dict(_PARAMS)
    for key, value in raw.items():
        mapped = _ALIASES.get(key, key)
        if mapped in params:
            params[mapped] = value
        elif key not in ("data", "counts", "subsample_indices", "seed_rng"):
            raise ValueError(f"unknown parameter: {key!r}")
    return params


def run(payload: dict) -> dict:
    """Run SPRING on a JSON payload and return a JSON-serialisable result."""
    if not isinstance(payload, dict):
        raise ValueError("input JSON must be an object.")
    data = payload.get("data", payload.get("counts"))
    if data is None:
        raise ValueError("input JSON must contain a 'data' (or 'counts') matrix.")
    matrix = validate_counts(data)

    params = _normalise_params(payload)
    subsample_indices = payload.get("subsample_indices")
    if subsample_indices is not None:
        subsample_indices = [np.asarray(ix, dtype=int) for ix in subsample_indices]

    result = spring(matrix, subsample_indices=subsample_indices, **params)

    stars = result.output.stars
    opt = stars.opt_index
    adjacency = np.asarray(result.adjacency, dtype=float)
    beta = np.asarray(result.beta.beta[opt - 1], dtype=float)
    weights = (beta + beta.T) / 2.0 * adjacency

    return {
        "spring_version": __version__,
        "lambdaseq": result.lambdaseq.tolist(),
        "opt_index": int(opt),
        "lambda": float(result.lambdaseq[opt - 1]),
        "instability": stars.summary.tolist(),
        "selection_frequency": np.asarray(stars.merge[opt - 1]).tolist(),
        "adjacency": adjacency.tolist(),
        "beta": beta.tolist(),
        "weights": weights.tolist(),
        "n_edges": int(adjacency.sum() / 2),
        "n": int(matrix.shape[0]),
        "p": int(matrix.shape[1]),
        "params": {
            k: (v if not isinstance(v, np.generic) else v.item())
            for k, v in params.items()
            if k != "lambdaseq"
        },
    }


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="spring-microbiome",
        description="SPRING network estimation. Reads JSON, writes JSON.",
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="-",
        help="input JSON file, or '-' for stdin (default)",
    )
    parser.add_argument(
        "-o", "--output", default="-", help="output JSON file, or '-' for stdout"
    )
    parser.add_argument("--pretty", action="store_true", help="pretty-print JSON")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args(argv)

    try:
        if args.input == "-":
            payload = json.load(sys.stdin)
        else:
            with open(args.input, encoding="utf-8") as fh:
                payload = json.load(fh)
        out = run(payload)
    except Exception as exc:  # noqa: BLE001 - CLI boundary
        print(f"error: {exc}", file=sys.stderr)
        return 1

    text = json.dumps(out, indent=2 if args.pretty else None)
    if args.output == "-":
        sys.stdout.write(text + "\n")
    else:
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(text + "\n")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
