"""Synthetic data generation for tests and validation."""

from __future__ import annotations

import numpy as np
import pandas as pd


def chain_precision_matrix(p: int, rho: float = 0.25) -> np.ndarray:
    """Create a sparse positive-definite chain precision matrix."""
    precision = np.eye(p)
    for i in range(p - 1):
        precision[i, i + 1] = -rho
        precision[i + 1, i] = -rho
    precision += 0.2 * np.eye(p)
    return precision


def _force_diag_dominance(precision: np.ndarray, margin: float = 0.25) -> np.ndarray:
    """Enforce strict diagonal dominance to keep the matrix positive definite."""
    A = precision.copy()
    abs_row_sum = np.sum(np.abs(A), axis=1) - np.abs(np.diag(A))
    np.fill_diagonal(A, abs_row_sum + margin)
    return A


def block_precision_matrix(
    p: int,
    *,
    n_blocks: int = 3,
    rho_intra: float = 0.22,
    rho_inter: float = 0.04,
) -> np.ndarray:
    """Create a block-structured sparse precision matrix."""
    if p < n_blocks:
        raise ValueError(f"p must be >= n_blocks (got p={p}, n_blocks={n_blocks}).")

    precision = np.zeros((p, p), dtype=float)
    block_size = int(np.floor(p / n_blocks))

    starts: list[int] = []
    ends: list[int] = []
    start = 0
    for b in range(n_blocks):
        extra = 1 if b < (p % n_blocks) else 0
        end = start + block_size + extra
        starts.append(start)
        ends.append(end)
        start = end

    for b in range(n_blocks):
        i0, i1 = starts[b], ends[b]
        for i in range(i0, i1):
            for j in range(i + 1, i1):
                precision[i, j] = -rho_intra
                precision[j, i] = -rho_intra

    for b in range(n_blocks - 1):
        left = ends[b] - 1
        right = starts[b + 1]
        precision[left, right] = -rho_inter
        precision[right, left] = -rho_inter

    return _force_diag_dominance(precision)


def hub_precision_matrix(
    p: int,
    *,
    n_hubs: int = 2,
    rho_hub: float = 0.20,
    rho_backbone: float = 0.06,
) -> np.ndarray:
    """Create a sparse hub-and-spoke precision matrix."""
    if p < 4:
        raise ValueError(f"p must be >= 4 (got p={p}).")
    if n_hubs < 1:
        raise ValueError(f"n_hubs must be >= 1 (got {n_hubs}).")

    precision = np.zeros((p, p), dtype=float)
    hubs = np.linspace(0, p - 1, num=n_hubs, dtype=int)

    for i in range(p):
        closest = hubs[int(np.argmin(np.abs(hubs - i)))]
        if closest != i:
            precision[i, closest] = -rho_hub
            precision[closest, i] = -rho_hub

    for i in range(p - 1):
        precision[i, i + 1] = min(precision[i, i + 1], -rho_backbone)
        precision[i + 1, i] = min(precision[i + 1, i], -rho_backbone)

    return _force_diag_dominance(precision)


def adjacency_from_precision(precision: np.ndarray) -> np.ndarray:
    """Create undirected adjacency from a precision matrix."""
    A = (np.abs(precision) > 1e-12).astype(int)
    np.fill_diagonal(A, 0)
    return A


def generate_counts_from_precision(
    precision: np.ndarray,
    *,
    n_samples: int,
    depth: int,
    seed: int,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Sample compositional counts from a precision matrix."""
    rng = np.random.default_rng(seed)
    covariance = np.linalg.inv(precision)
    n_taxa = int(precision.shape[0])

    latent = rng.multivariate_normal(np.zeros(n_taxa), covariance, size=n_samples)
    exp_latent = np.exp(latent)
    probs = exp_latent / np.sum(exp_latent, axis=1, keepdims=True)

    counts = np.vstack([rng.multinomial(depth, p) for p in probs])
    columns = [f"taxon_{i}" for i in range(n_taxa)]
    df = pd.DataFrame(counts, columns=columns)
    return df, adjacency_from_precision(precision)


def generate_compositional_counts(
    *,
    n_samples: int = 80,
    n_taxa: int = 8,
    depth: int = 2_000,
    rho: float = 0.25,
    seed: int = 7,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Generate synthetic count table plus ground-truth adjacency."""
    precision = chain_precision_matrix(n_taxa, rho=rho)
    return generate_counts_from_precision(
        precision,
        n_samples=n_samples,
        depth=depth,
        seed=seed,
    )


def generate_strong_signal_counts(
    *,
    n_samples: int = 220,
    n_taxa: int = 24,
    depth: int = 6_000,
    rho: float = 0.45,
    seed: int = 7,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Generate higher-signal compositional counts for robust parity testing."""
    return generate_compositional_counts(
        n_samples=n_samples,
        n_taxa=n_taxa,
        depth=depth,
        rho=rho,
        seed=seed,
    )


def generate_block_signal_counts(
    *,
    n_samples: int = 220,
    n_taxa: int = 24,
    depth: int = 6_000,
    n_blocks: int = 3,
    rho_intra: float = 0.24,
    rho_inter: float = 0.05,
    seed: int = 7,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Generate block-structured synthetic compositional counts."""
    precision = block_precision_matrix(
        n_taxa,
        n_blocks=n_blocks,
        rho_intra=rho_intra,
        rho_inter=rho_inter,
    )
    return generate_counts_from_precision(
        precision,
        n_samples=n_samples,
        depth=depth,
        seed=seed,
    )


def generate_hub_signal_counts(
    *,
    n_samples: int = 220,
    n_taxa: int = 24,
    depth: int = 6_000,
    n_hubs: int = 2,
    rho_hub: float = 0.20,
    rho_backbone: float = 0.06,
    seed: int = 7,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Generate hub-structured synthetic compositional counts."""
    precision = hub_precision_matrix(
        n_taxa,
        n_hubs=n_hubs,
        rho_hub=rho_hub,
        rho_backbone=rho_backbone,
    )
    return generate_counts_from_precision(
        precision,
        n_samples=n_samples,
        depth=depth,
        seed=seed,
    )
