"""Command-line interface: JSON in (stdin or file) -> JSON out (stdout or file).

Contract (see docs/PACKAGING.md):
- reads one JSON object from stdin or ``--input``;
- writes one JSON object to stdout or ``--output``; **nothing else** on stdout;
- warnings and errors go to stderr; errors exit with a non-zero status.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any

import numpy as np

from ._core import cclasso

__all__ = ["main", "run_request"]

# Request keys that map straight to cclasso() parameters.
_PARAM_KEYS = (
    "counts",
    "pseudo",
    "k_cv",
    "lam_int",
    "k_max",
    "n_boot",
    "seed",
)


def _array_or_none(req: dict, key: str, dtype) -> np.ndarray | None:
    if key not in req or req[key] is None:
        return None
    return np.asarray(req[key], dtype=dtype)


def run_request(req: dict[str, Any]) -> dict[str, Any]:
    """Execute a single request dict and return a JSON-serialisable response dict."""
    if "x" not in req:
        raise ValueError("request is missing required key 'x'.")
    x = np.asarray(req["x"], dtype=float)

    kwargs = {k: req[k] for k in _PARAM_KEYS if k in req and req[k] is not None}
    if "lam_int" in kwargs:
        kwargs["lam_int"] = tuple(kwargs["lam_int"])

    boot_indices = _array_or_none(req, "boot_indices", np.intp)
    null_norm = _array_or_none(req, "null_norm", float)

    res = cclasso(
        x,
        boot_indices=boot_indices,
        null_norm=null_norm,
        **kwargs,
    )
    return _to_jsonable(res)


def _to_jsonable(res: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key in ("var_w", "cor_w", "p_vals", "lambda"):
        out[key] = _encode(res[key])
    info = res["info_cv"]
    out["info_cv"] = {
        "lams": _encode(info["lams"]),
        "fvals": _encode(info["fvals"]),
        "k": int(info["k"]),
        "lam_int": _encode(info["lam_int"]),
    }
    return out


def _encode(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="cclasso",
        description="Python port of CCLasso (compositional correlation via Lasso).",
    )
    parser.add_argument("--input", "-i", help="input JSON file (default: stdin)")
    parser.add_argument("--output", "-o", help="output JSON file (default: stdout)")
    parser.add_argument("--indent", type=int, default=None, help="JSON indent for output")
    args = parser.parse_args(argv)

    try:
        if args.input:
            with open(args.input, encoding="utf-8") as fh:
                req = json.load(fh)
        else:
            req = json.load(sys.stdin)

        response = run_request(req)

        text = json.dumps(response, indent=args.indent)
        if args.output:
            with open(args.output, "w", encoding="utf-8") as fh:
                fh.write(text)
        else:
            sys.stdout.write(text)
            sys.stdout.write("\n")
    except Exception as exc:  # noqa: BLE001 - CLI boundary: report and exit non-zero
        sys.stderr.write(f"cclasso: error: {exc}\n")
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
