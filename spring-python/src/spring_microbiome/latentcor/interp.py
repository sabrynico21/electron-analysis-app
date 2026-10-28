"""Multilinear interpolation on a non-uniform Cartesian grid.

Faithful port of ``latentcor/src/ipol.c`` (``C_evalmlip``, linear blend),
including the boundary conventions of ``binsearch`` (never index 0, clamp to
``n-1``) and linear *extrapolation* outside the box.

The interpolation tables shipped in ``data/ipol_22.npz`` were extracted from
``latentcor`` 2.0.2 (see ``scripts/extract_ipol_tables.R``).
"""

from __future__ import annotations

from functools import lru_cache
from importlib import resources

import numpy as np

__all__ = ["MultilinearInterpolant", "load_ipol22"]


class MultilinearInterpolant:
    """Multilinear interpolant on a Cartesian product of 1-D grids.

    Parameters
    ----------
    grids : list of 1-D numpy.ndarray
        Ordered (strictly increasing) grid points per dimension.
    values : numpy.ndarray
        Array with ``values.shape == tuple(len(g) for g in grids)``.
    """

    def __init__(self, grids, values):
        self.grids = [np.asarray(g, dtype=float) for g in grids]
        self.values = np.asarray(values, dtype=float)
        expected = tuple(g.size for g in self.grids)
        if self.values.shape != expected:
            raise ValueError(f"values shape {self.values.shape} != grids {expected}")
        self.rank = len(self.grids)

    def _locate(self, dim: int, x: np.ndarray):
        g = self.grids[dim]
        n = g.size
        gp = np.searchsorted(g, x, side="left")
        gp = np.clip(gp, 1, n - 1)
        w = (x - g[gp - 1]) / (g[gp] - g[gp - 1])
        return gp, w

    def __call__(self, x: np.ndarray) -> np.ndarray:
        """Evaluate at columns of ``x`` (shape ``(rank, m)``)."""
        x = np.atleast_2d(np.asarray(x, dtype=float))
        if x.shape[0] != self.rank:
            x = x.T
        if x.shape[0] != self.rank:
            raise ValueError(f"expected {self.rank} coordinates per point")

        gp = np.empty((self.rank, x.shape[1]), dtype=np.intp)
        w = np.empty((self.rank, x.shape[1]), dtype=float)
        for d in range(self.rank):
            gp[d], w[d] = self._locate(d, x[d])

        out = np.zeros(x.shape[1], dtype=float)
        wsum = np.zeros(x.shape[1], dtype=float)
        n_corners = 1 << self.rank
        for corner in range(n_corners):
            idx = []
            cw = np.ones(x.shape[1], dtype=float)
            for d in range(self.rank):
                if (1 << d) & corner:
                    idx.append(gp[d])
                    cw = cw * w[d]
                else:
                    idx.append(gp[d] - 1)
                    cw = cw * (1.0 - w[d])
            out = out + cw * self.values[tuple(idx)]
            wsum = wsum + cw
        return out / wsum


@lru_cache(maxsize=1)
def load_ipol22() -> MultilinearInterpolant:
    """Load the shipped ``ipol_22`` table (Kendall x zratio x zratio)."""
    with resources.files("spring_microbiome.latentcor").joinpath(
        "data/ipol_22.npz"
    ).open("rb") as fh:
        with np.load(fh) as data:
            grid_k = data["grid_k"]
            grid_z1 = data["grid_z1"]
            grid_z2 = data["grid_z2"]
            values = data["values"]
            scale = float(data["scale"])
    interp = MultilinearInterpolant([grid_k, grid_z1, grid_z2], values)
    interp.scale = scale  # type: ignore[attr-defined]
    return interp
