#!/usr/bin/env python3
"""Write machine/environment metadata for an equivalence-fixture directory."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from pathlib import Path


def _git_revision(root: Path) -> dict[str, object]:
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "diff", "--quiet"], cwd=root, check=False
        ).returncode != 0
        return {"revision": revision, "tracked_worktree_dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"revision": None, "tracked_worktree_dirty": None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--component", required=True)
    args = parser.parse_args()

    try:
        import numpy as np
    except ImportError as exc:
        raise SystemExit(f"NumPy is required to write the manifest: {exc}") from exc

    packages = (
        "numpy",
        "scipy",
        "scikit-learn",
        "networkx",
        "pytest",
        "threadpoolctl",
    )
    versions = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None

    try:
        from threadpoolctl import threadpool_info

        threadpools = threadpool_info()
    except ImportError:
        threadpools = None

    numpy_config = np.__config__.show(mode="dicts")
    build_deps = numpy_config.get("Build Dependencies", {})
    blas_lapack = {
        key: {
            field: value.get(field)
            for field in ("name", "version", "openblas configuration", "found")
        }
        for key, value in build_deps.items()
        if key in {"blas", "lapack"} and isinstance(value, dict)
    }

    root = Path(__file__).resolve().parents[1]
    manifest = {
        "component": args.component,
        "python": {
            "version": sys.version,
            "implementation": platform.python_implementation(),
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "platform": platform.platform(),
        },
        "packages": versions,
        "thread_environment": {
            name: os.environ.get(name)
            for name in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            )
        },
        "threadpools": (
            [
                {
                    key: pool.get(key)
                    for key in ("user_api", "internal_api", "version", "num_threads", "architecture")
                }
                for pool in threadpools
            ]
            if threadpools is not None
            else None
        ),
        "blas_lapack": blas_lapack,
        "git": _git_revision(root),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
