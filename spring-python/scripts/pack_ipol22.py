#!/usr/bin/env python3
"""Pack the raw little-endian dumps produced by ``extract_ipol_tables.R``.

The result is a single ``ipol_22.npz`` shipped with the package.

Run:

    Rscript scripts/extract_ipol_tables.R /tmp/ipol_raw
    python scripts/pack_ipol22.py /tmp/ipol_raw src/spring_microbiome/latentcor/data/ipol_22.npz
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> None:
    """Read the raw dumps and write the packed ``.npz`` file."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw_dir", type=Path, help="directory with the .f64/.i32 dumps")
    parser.add_argument("out", type=Path, help="output .npz path")
    args = parser.parse_args()

    raw = args.raw_dir
    meta = json.loads((raw / "ipol_meta.json").read_text())
    dims = tuple(int(d) for d in meta["dims"])
    scale = float(meta["scale"])

    grid_k = np.fromfile(raw / "ipol_22_grid_K.f64", dtype="<f8")
    grid_z1 = np.fromfile(raw / "ipol_22_grid_z1.f64", dtype="<f8")
    grid_z2 = np.fromfile(raw / "ipol_22_grid_z2.f64", dtype="<f8")
    values_flat = np.fromfile(raw / "ipol_22_values.i32", dtype="<i4")

    assert grid_k.size == dims[0], (grid_k.size, dims)
    assert grid_z1.size == dims[1], (grid_z1.size, dims)
    assert grid_z2.size == dims[2], (grid_z2.size, dims)
    assert values_flat.size == int(np.prod(dims)), (values_flat.size, dims)

    # R arrays are column-major: first dimension varies fastest.
    values = values_flat.reshape(dims, order="F").astype(np.float64)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        grid_k=grid_k,
        grid_z1=grid_z1,
        grid_z2=grid_z2,
        values=values,
        scale=np.float64(scale),
    )
    print(f"wrote {args.out} ({args.out.stat().st_size} bytes), dims={dims}, scale={scale:g}")


if __name__ == "__main__":
    main()
