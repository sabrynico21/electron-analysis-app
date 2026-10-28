#!/usr/bin/env python3
"""
Main analysis script.
Reads parameters from --params JSON file, performs the analysis,
writes results to a JSON file in the OS temp directory.

Progress is reported to stdout as JSON lines: {"progress": <0-100>}
Regular log output goes to stdout as plain text.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

# Keep the statistical ports as standalone source packages and import them from
# their development/package locations. PyInstaller builds include these paths
# and their hidden imports in backend.spec.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _PORT_SOURCE in (
    os.path.join(_REPO_ROOT, 'cclasso-python', 'src'),
    os.path.join(_REPO_ROOT, 'spiec-easi-python', 'src'),
    os.path.join(_REPO_ROOT, 'spring-python', 'src'),
):
    if os.path.isdir(_PORT_SOURCE) and _PORT_SOURCE not in sys.path:
        sys.path.insert(0, _PORT_SOURCE)

P_VALUE_METHODS = {'spearman', 'sparcc', 'cclasso'}
NETWORK_METHODS = {'spieceasi', 'spring'}
CORRELATION_METHODS = P_VALUE_METHODS | NETWORK_METHODS
MATRIX_METHODS = {'sparcc', 'cclasso'} | NETWORK_METHODS

def report_progress(pct: int):
    print(json.dumps({"progress": pct}), flush=True)

def log(msg: str):
    print(msg, flush=True)


def _slug(value: str) -> str:
    return ''.join(ch.lower() if ch.isalnum() else '_' for ch in value).strip('_') or 'condition'


def _method_settings(method: str, params: dict) -> dict:
    """Return the effective model settings for result provenance."""
    if method == 'cclasso':
        # CCLasso is fed the same abundance table as every other method, so the
        # count -> compositional closure (0.5 pseudo-count) is always applied.
        return {
            'counts': True,
            'pseudo_count': float(params.get('cclassoPseudo', 0.5)),
            'k_cv': int(params.get('cclassoFolds', 3)),
            'n_boot': int(params.get('cclassoBootstraps', 20)),
            'seed': int(params.get('cclassoSeed', 17)),
        }
    if method == 'spieceasi':
        return {
            'method': 'mb',
            'selection': 'stars',
            'nlambda': int(params.get('spieceasiNLambdas', 20) or 20),
            'lambda_min_ratio': float(params.get('spieceasiLambdaMinRatio', 1e-3) or 1e-3),
            'rep_num': int(params.get('spieceasiNSubsamples', 20) or 20),
            'stars_threshold': float(params.get('spieceasiStarsThreshold', 0.05) or 0.05),
            'subsample_fraction': params.get('spieceasiSubsampleFraction'),
            'seed': 17,
        }
    if method == 'spring':
        # SPRING receives the same abundance table as every other method and
        # applies its own compositional transform (MCLR) internally.
        return {
            'method': 'mb',
            'quantitative': False,
            'lambda_path': 'data-specific',
            'lambda_min_ratio': float(params.get('springLambdaMinRatio', 1e-2)),
            'nlambda': int(params.get('springNLambdas', 20) or 20),
            'rep_num': int(params.get('springNSubsamples', 20) or 20),
            'stars_threshold': float(params.get('springStarsThreshold', 0.1)),
            'subsample_fraction': float(params.get('springSubsampleFraction', 0.8)),
            'seed': int(params.get('springSeed', 17)),
            'Rmethod': str(params.get('springRmethod', 'approx')),
            'use_nearPD': bool(params.get('springUseNearPD', True)),
            'nu': float(params.get('springNu', 0.001)),
            'ratio': float(params.get('springRatio', 0.9)),
        }
    if method == 'sparcc':
        return {'permutations': int(params.get('sparccPermutations', 100) or 100)}
    return {'method': 'spearman'}


def _detect_delimiter(header_line: str) -> str:
    if '\t' in header_line:
        return '\t'
    if ';' in header_line:
        return ';'
    return ','


def _read_table_rows(file_path: str) -> list:
    # Windows defaults to the ANSI codepage, and Excel exports usually carry a BOM:
    # utf-8-sig handles both without mangling the first column name. `errors`
    # degrades gracefully instead of aborting the whole analysis on one bad byte.
    with open(file_path, 'r', newline='', encoding='utf-8-sig', errors='replace') as handle:
        first_line = handle.readline()
        if not first_line:
            return []
        delimiter = _detect_delimiter(first_line)
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter=delimiter)
        return [row for row in reader]


def _load_abundance_matrices(file_paths: list, known_sample_ids: set[str] | None = None) -> tuple[dict, list]:
    """
    Returns:
      sample_to_values: { sample_id: { feature_name: numeric_value } }
      ordered_features: unique feature names in first-seen order
    """
    sample_to_values = {}
    ordered_features = []
    seen_features = set()

    for file_path in file_paths:
        rows = _read_table_rows(file_path)
        if not rows:
            continue

        keys = list(rows[0].keys())
        id_key = keys[0]
        other_keys = [key for key in keys if key != id_key]

        # Detect orientation. Default is samples on rows, features on columns.
        samples_on_columns = False
        if known_sample_ids:
            row_id_matches = 0
            for row in rows:
                row_id = (row.get(id_key) or '').strip()
                if row_id in known_sample_ids:
                    row_id_matches += 1

            col_id_matches = 0
            for key in other_keys:
                if (key or '').strip() in known_sample_ids:
                    col_id_matches += 1

            # If headers look like sample IDs more than first-column values, table is transposed.
            samples_on_columns = col_id_matches > row_id_matches and col_id_matches > 0

        if samples_on_columns:
            log(f"Detected transposed abundance matrix in '{file_path}': samples on columns, features on rows")
            for row in rows:
                feature_name = (row.get(id_key) or '').strip()
                if not feature_name:
                    continue
                if feature_name not in seen_features:
                    seen_features.add(feature_name)
                    ordered_features.append(feature_name)

                for sample_id in other_keys:
                    sample_name = (sample_id or '').strip()
                    if not sample_name:
                        continue
                    raw_value = (row.get(sample_id) or '').strip()
                    if raw_value == '':
                        continue
                    try:
                        numeric_value = float(raw_value)
                    except ValueError:
                        continue
                    sample_values = sample_to_values.setdefault(sample_name, {})
                    sample_values[feature_name] = numeric_value
        else:
            sample_key = id_key
            feature_names = other_keys

            for feature_name in feature_names:
                if feature_name not in seen_features:
                    seen_features.add(feature_name)
                    ordered_features.append(feature_name)

            for row in rows:
                sample_id = (row.get(sample_key) or '').strip()
                if not sample_id:
                    continue
                sample_values = sample_to_values.setdefault(sample_id, {})
                for feature_name in feature_names:
                    raw_value = (row.get(feature_name) or '').strip()
                    if raw_value == '':
                        continue
                    try:
                        sample_values[feature_name] = float(raw_value)
                    except ValueError:
                        continue

    return sample_to_values, ordered_features


def _load_metadata_conditions(
    metadata_paths: list,
    domain_label: str = 'general',
    condition_column: str = 'Combination_treat2',
) -> dict:
    """Returns {sample_id: condition_value} from metadata files.

    Args:
        metadata_paths: List of metadata file paths
        domain_label: Label for logging (e.g., 'bacteria', 'fungi', 'general')
        condition_column: Name of the column holding the experimental condition
    """
    sample_to_condition = {}
    for file_path in metadata_paths:
        rows = _read_table_rows(file_path)
        if not rows:
            continue

        keys = list(rows[0].keys())
        if condition_column not in keys:
            raise ValueError(f"Metadata file '{file_path}' does not contain {condition_column}")
        sample_key = keys[0]

        for row in rows:
            sample_id = (row.get(sample_key) or '').strip()
            condition = (row.get(condition_column) or '').strip()
            if sample_id and condition:
                sample_to_condition[sample_id] = condition

    log(f"  {domain_label} metadata: {len(sample_to_condition)} samples loaded")
    return sample_to_condition


def _compute_spearman(x_values: list, y_values: list) -> tuple[float, float]:
    # Imported here rather than at module import: scipy.stats is heavy (about a
    # second in the frozen bundle) and this helper is only reached by the
    # Pearson/Spearman path, never by the SparCC or SPIEC-EASI runs.
    from scipy.stats import spearmanr

    # Spearman is undefined for constant vectors and may emit warnings; treat as invalid pair.
    if len(x_values) < 3 or len(y_values) < 3:
        return math.nan, math.nan

    if len(set(x_values)) <= 1 or len(set(y_values)) <= 1:
        return math.nan, math.nan

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        corr, p_value = spearmanr(x_values, y_values)

    return float(corr), float(p_value)


def _record_method_drop(
    drop_events: list | None,
    method: str,
    scope: str,
    kind: str,
    items: list,
    samples_in_scope: int,
) -> None:
    """Record a removal performed by a correlation method itself.

    The user filters (prevalence / mean relative abundance / variability) run
    first and are fully reported; the backends then still discard taxa that are
    constant inside the analysed sample subset (and SPRING discards all-zero
    samples). That removal used to exist only as a log line, which made the
    reported ``before``/``after`` counts disagree with the number of taxa the
    model actually saw. Events collected here are attached to the
    ``taxa_filtering`` report of the results payload.

    ``kind`` is ``'constant_feature'`` (items are taxa) or ``'empty_sample'``
    (items are sample ids); ``scope`` names the condition/matrix that was fitted
    (e.g. ``'bacteria · control'``, ``'bac+fun · treated'``).
    """
    if drop_events is None or not items:
        return
    preview = list(items[:50])
    drop_events.append({
        'method': method,
        'scope': scope,
        'kind': kind,
        'item_kind': 'taxa' if kind == 'constant_feature' else 'samples',
        'count': len(items),
        'items': preview,
        'items_truncated': len(items) > len(preview),
        'samples_in_scope': int(samples_in_scope),
    })


def _matrix_for_samples(
    matrix_by_sample: dict,
    features: list,
    sample_ids: list,
    drop_events: list | None = None,
    method: str = '',
    scope: str = '',
) -> tuple[list[str], np.ndarray]:
    """Return a finite samples × features matrix with constant taxa removed."""
    usable_samples = [sample_id for sample_id in sample_ids if sample_id in matrix_by_sample]
    varying_features, dropped = _filter_constant_features(
        features, usable_samples, matrix_by_sample
    )
    if dropped:
        log(
            f"  Removed {len(dropped)} constant feature(s) before model fitting "
            f"(zero variance / all-zero across samples)."
        )
        _record_method_drop(
            drop_events, method, scope, 'constant_feature', dropped, len(usable_samples)
        )
    matrix = np.asarray(
        [
            [float(matrix_by_sample[sample_id].get(feature, 0.0)) for feature in varying_features]
            for sample_id in usable_samples
        ],
        dtype=float,
    )
    if matrix.size and not np.isfinite(matrix).all():
        raise ValueError('The selected abundance matrix contains non-finite values.')
    return varying_features, matrix


def _compute_cclasso_for_matrix(
    matrix_by_sample: dict,
    features: list,
    sample_ids: list,
    params: dict,
    progress_range: tuple[float, float] | None = None,
    drop_events: list | None = None,
    scope_label: str = '',
) -> tuple[dict, dict, dict]:
    """Infer composition-aware correlations and bootstrap p-values with CCLasso."""
    from cclasso import cclasso

    start_pct, end_pct = (progress_range or (0.0, 0.0))

    def _progress(fraction: float) -> None:
        if progress_range is not None:
            report_progress(round(start_pct + (end_pct - start_pct) * min(1.0, max(0.0, fraction))))

    usable_samples = [sample_id for sample_id in sample_ids if sample_id in matrix_by_sample]
    varying_features, X = _matrix_for_samples(
        matrix_by_sample, features, usable_samples, drop_events, 'cclasso', scope_label
    )
    if len(usable_samples) < 6 or len(varying_features) < 2:
        raise ValueError(
            'CCLasso needs at least 6 samples and 2 non-constant taxa '
            f'(got n={len(usable_samples)}, p={len(varying_features)}).'
        )
    if np.any(X < 0):
        raise ValueError('CCLasso requires non-negative abundance values.')

    # Same input as every other correlation method: the abundance table as loaded.
    # CCLasso's own model turns it into compositions (pseudo-count + row closure).
    counts = True
    pseudo = float(params.get('cclassoPseudo', 0.5))
    k_cv = int(params.get('cclassoFolds', 3))
    n_boot = int(params.get('cclassoBootstraps', 20))
    seed = int(params.get('cclassoSeed', 17))

    _progress(0.05)
    log(
        f"Running CCLasso: {len(usable_samples)} samples × {len(varying_features)} taxa "
        f"(counts={counts}, k_cv={k_cv}, n_boot={n_boot})"
    )
    try:
        with warnings.catch_warnings(record=True) as solver_warnings:
            warnings.simplefilter('always')
            result = cclasso(
                X,
                counts=counts,
                pseudo=pseudo,
                k_cv=k_cv,
                n_boot=n_boot,
                seed=seed,
            )
    except Exception as exc:
        raise RuntimeError(f'CCLasso correlation inference failed: {exc}') from exc

    # Solver warnings (ADMM iteration cap, lambda near boundary, numpy numerics)
    # are captured instead of being dumped to stderr as raw RuntimeWarning text,
    # which the log viewer would misclassify. They stay reported, as warnings.
    for item in solver_warnings:
        message = ' '.join(str(item.message).split())
        log(f'  CCLasso solver note: {message}')
    _progress(1.0)

    correlations = result['cor_w']
    p_values = result['p_vals']
    corr_matrix = {feature: {} for feature in varying_features}
    pval_matrix = {feature: {} for feature in varying_features}
    finite_pairs = 0
    for i, left in enumerate(varying_features):
        for j in range(i + 1, len(varying_features)):
            right = varying_features[j]
            corr = float(correlations[i, j])
            p_value = float(p_values[i, j])
            if not (math.isfinite(corr) and math.isfinite(p_value)):
                continue
            corr_matrix[left][right] = corr
            corr_matrix[right][left] = corr
            pval_matrix[left][right] = p_value
            pval_matrix[right][left] = p_value
            finite_pairs += 1

    diagnostics = {
        'weight_type': 'correlation',
        'evidence_type': 'bootstrap_p_value',
        'lambda': float(result['lambda']),
        'n_samples': len(usable_samples),
        'n_taxa': len(varying_features),
        'finite_pair_count': finite_pairs,
        'nonfinite_pair_count': (len(varying_features) * (len(varying_features) - 1) // 2) - finite_pairs,
        'admm_iteration_cap_count': int(
            np.count_nonzero(np.asarray(result.get('_admm_iters', [])) >= 200)
        ),
    }
    log(
        f"  CCLasso selected lambda={diagnostics['lambda']:.6g}; "
        f"finite pair estimates={finite_pairs}; "
        f"ADMM calls at 200-iteration cap={diagnostics['admm_iteration_cap_count']}"
    )
    if diagnostics['nonfinite_pair_count']:
        log(
            f"  WARNING: CCLasso returned non-finite estimates for "
            f"{diagnostics['nonfinite_pair_count']} pair(s); those pairs are omitted."
        )
    return corr_matrix, pval_matrix, diagnostics


def _write_counts_tsv(path: str, samples: list[str], features: list[str], matrix_by_sample: dict) -> None:
    # FastSpar requires the classic BIOM TSV format: the first header cell must be
    # '#OTU ID', each row is an OTU (feature) and each column is a sample.
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['#OTU ID'] + samples)
        for feature in features:
            row = [feature]
            for sample_id in samples:
                values = matrix_by_sample.get(sample_id, {})
                value = values.get(feature, 0.0)
                row.append(float(value) if value is not None else 0.0)
            writer.writerow(row)


def _read_square_matrix_tsv(path: str) -> tuple[list[str], dict]:
    with open(path, 'r', newline='', encoding='utf-8') as handle:
        reader = csv.reader(handle, delimiter='\t')
        rows = [row for row in reader if row]

    if not rows:
        return [], {}

    header = rows[0]
    if len(header) < 2:
        return [], {}

    features = [str(item).strip() for item in header[1:] if str(item).strip()]
    matrix = {feature: {} for feature in features}

    for row in rows[1:]:
        if len(row) < 2:
            continue
        row_feature = str(row[0]).strip()
        if not row_feature:
            continue
        if row_feature not in matrix:
            matrix[row_feature] = {}
        for idx, col_feature in enumerate(features, start=1):
            if idx >= len(row):
                continue
            cell = str(row[idx]).strip()
            if not cell:
                continue
            try:
                matrix[row_feature][col_feature] = float(cell)
            except ValueError:
                continue

    return features, matrix


def _resolve_fastspar_binary(params: dict) -> str:
    candidate = (params.get('fastsparPath') or '').strip()
    if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate

    env_candidate = (os.environ.get('FASTSPAR_BIN') or '').strip()
    if env_candidate and os.path.isfile(env_candidate) and os.access(env_candidate, os.X_OK):
        return env_candidate

    found = shutil.which('fastspar')
    if found:
        return found

    raise RuntimeError(
        'SparCC requested but FastSpar binary is not available. '
        'Configure fastsparPath in Settings, set FASTSPAR_BIN, or install fastspar in PATH.'
    )


def _resolve_fastspar_sibling(fastspar_bin: str, name: str) -> str:
    """Resolve a FastSpar companion binary (e.g. fastspar_bootstrap or
    fastspar_pvalues) next to ``fastspar`` first, then fall back to PATH."""
    fastspar_dir = os.path.dirname(fastspar_bin)
    exe = name + ('.exe' if os.name == 'nt' else '')
    local = os.path.join(fastspar_dir, exe)
    if os.path.isfile(local) and os.access(local, os.X_OK):
        return local
    found = shutil.which(name)
    if found:
        return found
    return ''


def _run_and_stream(command: list, label: str, quiet: bool = False) -> tuple[int, str]:
    """Run a subprocess, forwarding its stdout+stderr lines into the analysis log.

    This makes FastSpar's own progress/output appear live in the "Run Analysis"
    log window. When ``quiet`` is True (used for the many bootstrap runs) only
    the failure tail is retained. Returns (returncode, tail_of_output) so
    failures can be reported with the last lines of the process output.
    """
    if not quiet:
        log(f"[{label}] $ {' '.join(command)}")
    start = time.perf_counter()
    proc = subprocess.Popen(
        command,
        stdin=subprocess.DEVNULL,  # never block on interactive stdin prompts
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding='utf-8',          # Windows would otherwise decode as ANSI
        errors='replace',
        bufsize=1,
    )
    tail_lines: list = []
    assert proc.stdout is not None
    for line in proc.stdout:
        line = line.rstrip('\n')
        if line.strip():
            if not quiet:
                log(f"[{label}] {line}")
            tail_lines.append(line)
    proc.wait()
    if not quiet:
        elapsed = time.perf_counter() - start
        log(f"[{label}] finished in {elapsed:.2f}s (exit code {proc.returncode})")
    return proc.returncode, '\n'.join(tail_lines[-20:])


def _filter_constant_features(
    features: list,
    usable_samples: list,
    matrix_by_sample: dict,
) -> tuple[list, list]:
    """Split features into (varying, constant) across the given samples.

    A feature whose value is identical in every sample (e.g. all-zero, or a
    fixed abundance) has only one unique permutation for FastSpar's p-value
    test and zero variance, so it cannot produce a meaningful correlation.
    These are exactly the OTUs FastSpar warns about; we drop them upstream.
    """
    varying = []
    constant = []
    for feature in features:
        values = [matrix_by_sample[sample].get(feature, 0.0) for sample in usable_samples]
        if len(set(values)) <= 1:
            constant.append(feature)
        else:
            varying.append(feature)
    return varying, constant


def _compute_sparcc_for_matrix(
    matrix_by_sample: dict,
    features: list,
    sample_ids: list,
    params: dict,
    progress_range: tuple[float, float] | None = None,
    drop_events: list | None = None,
    scope_label: str = '',
) -> tuple[dict, dict]:
    start_pct, end_pct = (progress_range or (0.0, 0.0))

    def _progress(fraction: float) -> None:
        if progress_range is not None:
            clamped = min(1.0, max(0.0, fraction))
            report_progress(round(start_pct + (end_pct - start_pct) * clamped))

    usable_samples = [sample_id for sample_id in sample_ids if sample_id in matrix_by_sample]
    if len(usable_samples) < 3 or len(features) < 2:
        return {}, {}
    _progress(0.0)

    # Drop constant features (same value in every sample) before handing the
    # matrix to FastSpar: they have zero variance, can't correlate with
    # anything, and FastSpar flags them as having "only one unique permutation".
    features, dropped = _filter_constant_features(features, usable_samples, matrix_by_sample)
    if dropped:
        preview = ', '.join(dropped[:20])
        if len(dropped) > 20:
            preview += f' … (+{len(dropped) - 20} more)'
        log(
            f"  Removed {len(dropped)} constant feature(s) before FastSpar "
            f"(zero variance / only one unique permutation): {preview}"
        )
    _record_method_drop(
        drop_events, 'sparcc', scope_label, 'constant_feature', dropped, len(usable_samples)
    )
    if len(features) < 2:
        log("  Fewer than 2 varying features remain; skipping SparCC matrix")
        return {}, {}

    fastspar_bin = _resolve_fastspar_binary(params)
    permutations = int(params.get('sparccPermutations', 100) or 100)
    permutations = max(10, permutations)

    # FastSpar's permutation p-values need the companion binaries: fastspar_bootstrap
    # (generates permuted count tables) and fastspar_pvalues (reads the bootstrap
    # correlations). Resolve them next to fastspar first, then fall back to PATH.
    pvalue_bin = _resolve_fastspar_sibling(fastspar_bin, 'fastspar_pvalues')
    bootstrap_bin = _resolve_fastspar_sibling(fastspar_bin, 'fastspar_bootstrap')
    missing = [name for name, path in (
        ('fastspar_pvalues', pvalue_bin),
        ('fastspar_bootstrap', bootstrap_bin),
    ) if not path]
    if missing:
        raise RuntimeError(
            'SparCC p-value computation requires the FastSpar companion binaries: '
            f'{", ".join(missing)}. Place them next to fastspar under '
            'resources/bin/<platform>/<arch>/.'
        )

    # Parallelism: default to one fastspar process per core (1 thread each) for
    # the bootstrap correlations. Both knobs can be overridden via params.
    threads = max(1, os.cpu_count() or 2)
    if isinstance(params.get('sparccThreads'), int) and params.get('sparccThreads', 0) > 0:
        threads = max(1, int(params['sparccThreads']))
    workers = min(permutations, max(1, threads))
    if isinstance(params.get('sparccWorkers'), int) and params.get('sparccWorkers', 0) > 0:
        workers = min(permutations, max(1, int(params['sparccWorkers'])))

    with tempfile.TemporaryDirectory(prefix='sparcc_') as tmp_dir:
        counts_path = os.path.join(tmp_dir, 'counts.tsv')
        median_corr_path = os.path.join(tmp_dir, 'median_correlation.tsv')
        median_cov_path = os.path.join(tmp_dir, 'median_covariance.tsv')
        bootstrap_prefix = os.path.join(tmp_dir, 'bootstrap_counts')
        corr_prefix = os.path.join(tmp_dir, 'bootstrap_corr')
        cov_prefix = os.path.join(tmp_dir, 'bootstrap_cov')
        pval_path = os.path.join(tmp_dir, 'pvalues.tsv')

        _write_counts_tsv(counts_path, usable_samples, features, matrix_by_sample)

        # 1) Median correlation on the original counts.
        log(
            f"Running SparCC correlation: {len(usable_samples)} samples x "
            f"{len(features)} features via FastSpar "
            f"({permutations} permutations for p-values)"
        )
        code, tail = _run_and_stream([
            fastspar_bin,
            '--otu_table', counts_path,
            '--correlation', median_corr_path,
            '--covariance', median_cov_path,
            '--iterations', '20',
            '--exclude_iterations', '10',
            '--threads', str(threads),
            '--yes',
        ], 'fastspar')
        if code != 0:
            raise RuntimeError(f'FastSpar failed (exit {code}): {tail or "unknown error"}')
        _progress(0.2)

        # 2) Generate <permutations> permuted count tables (the bootstrap).
        code, tail = _run_and_stream([
            bootstrap_bin,
            '--otu_table', counts_path,
            '--number', str(permutations),
            '--prefix', bootstrap_prefix,
            '--threads', str(threads),
        ], 'fastspar_bootstrap')
        if code != 0:
            raise RuntimeError(
                f'fastspar_bootstrap failed (exit {code}): {tail or "unknown error"}'
            )
        _progress(0.4)

        # 3) Infer a correlation matrix for every bootstrap counts table,
        #    running the independent runs in parallel (1 thread each).
        log(
            f"  Computing {permutations} bootstrap correlations in parallel "
            f"({workers} workers x 1 thread each; {threads} threads for single runs)…"
        )
        _bootstrap_start = time.perf_counter()

        def _bootstrap_correlation(i: int):
            return _run_and_stream([
                fastspar_bin,
                '--otu_table', f'{bootstrap_prefix}_{i}.tsv',
                '--correlation', f'{corr_prefix}_{i}.tsv',
                '--covariance', f'{cov_prefix}_{i}.tsv',
                '--iterations', '20',
                '--exclude_iterations', '10',
                '--threads', '1',
                '--yes',
            ], f'fastspar_bootstrap[{i}]', quiet=True)

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_bootstrap_correlation, i): i for i in range(permutations)}
            for done, future in enumerate(as_completed(futures), start=1):
                index = futures[future]
                code, tail = future.result()
                if code != 0:
                    raise RuntimeError(
                        f'fastspar bootstrap correlation {index} failed (exit {code}): '
                        f'{tail or "unknown error"}'
                    )
                elapsed = time.perf_counter() - _bootstrap_start
                log(
                    f"    SparCC bootstrap corr {done}/{permutations} "
                    f"({elapsed:.1f}s elapsed, {elapsed / done:.2f}s/run)"
                )
                _progress(0.4 + 0.55 * (done / permutations))

        log(
            f"  All {permutations} bootstrap correlations computed in "
            f"{time.perf_counter() - _bootstrap_start:.1f}s"
        )

        # 4) Compute permutation p-values from the bootstrap correlations.
        code, tail = _run_and_stream([
            pvalue_bin,
            '--otu_table', counts_path,
            '--correlation', median_corr_path,
            '--prefix', f'{corr_prefix}_',
            '--permutations', str(permutations),
            '--outfile', pval_path,
            '--threads', str(threads),
        ], 'fastspar_pvalues')
        if code != 0:
            raise RuntimeError(
                f'fastspar_pvalues failed (exit {code}): {tail or "unknown error"}'
            )
        _progress(1.0)

        _, pvalue_matrix = _read_square_matrix_tsv(pval_path)
        _, corr_matrix = _read_square_matrix_tsv(median_corr_path)
        return corr_matrix, pvalue_matrix


def _compute_spieceasi_for_matrix(
    matrix_by_sample: dict,
    features: list,
    sample_ids: list,
    params: dict,
    progress_range: tuple[float, float] | None = None,
    subsamples: list | None = None,
    drop_events: list | None = None,
    scope_label: str = '',
) -> tuple[dict, dict, dict]:
    """Compute SPIEC-EASI MB weights and expose them as correlation-like pairs.

    Edges selected by MB+StARS receive p=0.0 so they pass the p-value threshold;
    non-selected edges receive p=1.0 and are filtered out by default.

    ``subsamples`` optionally injects fixed 0-based StARS subsample index arrays,
    bypassing the internal RNG, so equivalence tests can feed R's exact
    subsamples. The default (``None``) keeps the normal random behaviour.
    """
    from spieceasi import spiec_easi

    start_pct, end_pct = (progress_range or (0.0, 0.0))

    def _progress(fraction: float) -> None:
        if progress_range is not None:
            clamped = min(1.0, max(0.0, fraction))
            report_progress(round(start_pct + (end_pct - start_pct) * clamped))

    usable_samples = [sample_id for sample_id in sample_ids if sample_id in matrix_by_sample]
    if len(usable_samples) < 4 or len(features) < 2:
        return {}, {}, {"selected_lambda": math.nan, "edge_count": 0}

    # Match the FastSpar path behavior: drop constant features (including
    # all-zero taxa) before model fitting. Such features carry no variance and
    # are not informative for network inference.
    features, dropped = _filter_constant_features(features, usable_samples, matrix_by_sample)
    if dropped:
        preview = ', '.join(dropped[:20])
        if len(dropped) > 20:
            preview += f' … (+{len(dropped) - 20} more)'
        log(
            f"  Removed {len(dropped)} constant feature(s) before SPIEC-EASI "
            f"(zero variance / all-zero across samples): {preview}"
        )
    _record_method_drop(
        drop_events, 'spieceasi', scope_label, 'constant_feature', dropped, len(usable_samples)
    )
    if len(features) < 2:
        log("  Fewer than 2 varying features remain; skipping SPIEC-EASI matrix")
        return {}, {}, {"selected_lambda": math.nan, "edge_count": 0}

    _progress(0.05)

    X = np.array(
        [
            [float(matrix_by_sample[sample_id].get(feature, 0.0)) for feature in features]
            for sample_id in usable_samples
        ],
        dtype=float,
    )
    n_samples, n_taxa = X.shape
    if n_samples <= n_taxa:
        log(
            "  WARNING: SPIEC-EASI MB with n_samples <= n_taxa can be unstable "
            f"(n_samples={n_samples}, n_taxa={n_taxa})."
        )

    # Defaults aligned with the R SpiecEasi reference (spiec.easi.default):
    # nlambda=20, lambda.min.ratio=1e-3, rep.num=20. Kept explicit so the app
    # and the R golden reference use exactly the same settings.
    n_lambdas = int(params.get('spieceasiNLambdas', 20) or 20)
    lambda_min_ratio = float(params.get('spieceasiLambdaMinRatio', 1e-3) or 1e-3)
    n_subsamples = int(params.get('spieceasiNSubsamples', 20) or 20)
    stars_threshold = float(params.get('spieceasiStarsThreshold', 0.05) or 0.05)
    # None -> follow pulsar::.ratcheck (ratio 0.8 for n <= 144, else 10*sqrt(n)/n),
    # matching the R SpiecEasi reference; an explicit value still overrides it.
    _raw_frac = params.get('spieceasiSubsampleFraction')
    subsample_fraction = float(_raw_frac) if _raw_frac else None
    # Worker processes for the StARS subsample loop. The port runs the
    # subsamples in a process pool and aggregates them in a fixed order, so
    # results are identical for any worker count. Default: ~physical cores.
    default_jobs = max(1, (os.cpu_count() or 1) // 2)
    n_jobs = int(params.get('spieceasiJobs', default_jobs) or 1)
    n_jobs = max(1, n_jobs)
    log(f"  SPIEC-EASI worker processes: {n_jobs}")

    # With very few samples and hundreds of taxa, MB+StARS can require a huge
    # number of Lasso fits. Reduce defaults in that regime unless the user has
    # explicitly overridden them.
    if n_samples <= 10 and n_taxa >= 200:
        # Trim the workload in the very-high-dimensional, tiny-sample regime
        # unless the user overrode the settings explicitly. (The port follows
        # R's fixed MB MAX_ITER=10000; there is no max_iter override.)
        if 'spieceasiNLambdas' not in params:
            n_lambdas = min(n_lambdas, 10)
        if 'spieceasiNSubsamples' not in params:
            n_subsamples = min(n_subsamples, 8)

    estimated_model_fits = (n_lambdas * n_subsamples + 1) * n_taxa
    log(
        "  SPIEC-EASI workload estimate: "
        f"{n_lambdas} lambdas x {n_subsamples} subsamples x {n_taxa} taxa "
        f"(+ final fit) ~= {estimated_model_fits:,} neighborhood Lasso fits"
    )
    log(
        "  SPIEC-EASI settings: "
        f"n_jobs={n_jobs}, "
        f"lambda_min_ratio={lambda_min_ratio}, stars_threshold={stars_threshold}"
    )

    def _on_progress(done: int, total: int) -> None:
        # Reserve 5% for setup and 95% for StARS+final fit.
        _progress(0.05 + 0.9 * (done / max(1, total)))
        log(f"  SPIEC-EASI StARS progress: {done}/{total} subsamples")

    try:
        result = spiec_easi(
            X,
            method='mb',
            sel_criterion='stars',
            verbose=False,
            nlambda=n_lambdas,
            lambda_min_ratio=lambda_min_ratio,
            pulsar_params={
                'rep_num': n_subsamples,
                'thresh': stars_threshold,
                'subsample_ratio': subsample_fraction,
                'seed': 17,
            },
            sym='or',
            subsample_indices=subsamples,
            n_jobs=n_jobs,
            progress_callback=_on_progress,
        )
    except Exception as exc:
        raise RuntimeError(f'SPIEC-EASI MB failed: {exc}') from exc

    _progress(1.0)

    adjacency = result['refit']
    edge_stability = result['merge'][result['opt_index'] - 1]
    selected_lambda = float(result['opt_lambda'])

    corr_matrix = {feature: {} for feature in features}
    pval_matrix = {feature: {} for feature in features}
    edge_count = 0
    for i, left in enumerate(features):
        for j in range(i + 1, len(features)):
            right = features[j]
            is_edge = int(adjacency[i, j]) != 0
            if not is_edge:
                continue
            # StARS selects network edges; it does not calculate pairwise
            # hypothesis-test p-values. Preserve its actual edge-frequency score
            # and leave p-value absent rather than fabricating 0/1 values.
            correlation = float(edge_stability[i, j])
            p_value = None
            edge_count += 1
            corr_matrix[left][right] = correlation
            corr_matrix[right][left] = correlation
            pval_matrix[left][right] = p_value
            pval_matrix[right][left] = p_value

    diagnostics = {
        'selected_lambda': selected_lambda,
        'edge_count': edge_count,
        'n_samples': int(X.shape[0]),
        'n_taxa': int(X.shape[1]),
        'mean_edge_stability': (
            float(np.mean(edge_stability[adjacency != 0])) if edge_count else 0.0
        ),
        'weight_type': 'selection_stability',
        'evidence_type': 'stars',
    }
    log(
        f"  SPIEC-EASI selected lambda={selected_lambda:.6g} with {edge_count} "
        f"StARS-selected edges"
    )
    return corr_matrix, pval_matrix, diagnostics


def _get_pair_stat(
    corr_matrix: dict,
    pval_matrix: dict,
    left: str,
    right: str,
) -> tuple[float, float | None]:
    corr = corr_matrix.get(left, {}).get(right)
    if corr is None:
        corr = corr_matrix.get(right, {}).get(left)
    pval = pval_matrix.get(left, {}).get(right)
    if pval is None:
        pval = pval_matrix.get(right, {}).get(left)

    if corr is None:
        return math.nan, None

    return float(corr), (None if pval is None else float(pval))


def _compute_spring_for_matrix(
    matrix_by_sample: dict,
    features: list,
    sample_ids: list,
    params: dict,
    progress_range: tuple[float, float] | None = None,
    drop_events: list | None = None,
    scope_label: str = '',
) -> tuple[dict, dict, dict]:
    """Run SPRING/StARS and return selected edges weighted by latent correlation.

    SPRING is a network selector, not a pairwise significance test. Only its
    StARS-selected edges are returned, their weights are the signed latent
    correlations, and their p-values remain ``None``.
    """
    from spring_microbiome import spring as spring_fit
    from spring_microbiome.latentcor import latentcor
    from spring_microbiome.mclr import mclr

    start_pct, end_pct = progress_range or (0.0, 0.0)

    def _progress(fraction: float) -> None:
        if progress_range is not None:
            report_progress(round(start_pct + (end_pct - start_pct) * min(1.0, max(0.0, fraction))))

    usable_samples = [sample_id for sample_id in sample_ids if sample_id in matrix_by_sample]
    varying_features, X = _matrix_for_samples(
        matrix_by_sample, features, usable_samples, drop_events, 'spring', scope_label
    )
    if np.any(X < 0):
        raise ValueError('SPRING requires non-negative abundance values.')
    nonzero_rows = X.sum(axis=1) > 0
    if not np.all(nonzero_rows):
        removed_samples = [
            sample for sample, keep in zip(usable_samples, nonzero_rows) if not keep
        ]
        log(f"  Removed {len(removed_samples)} all-zero sample(s) before SPRING.")
        _record_method_drop(
            drop_events, 'spring', scope_label, 'empty_sample',
            removed_samples, len(usable_samples),
        )
        X = X[nonzero_rows]
        usable_samples = [sample for sample, keep in zip(usable_samples, nonzero_rows) if keep]
        # Dropping samples can make further taxa constant: record that too.
        varying_features, X = _matrix_for_samples(
            matrix_by_sample, varying_features, usable_samples,
            drop_events, 'spring', scope_label,
        )
    if len(usable_samples) < 6 or len(varying_features) < 2:
        raise ValueError(
            'SPRING needs at least 6 samples and 2 non-constant taxa '
            f'(got n={len(usable_samples)}, p={len(varying_features)}).'
        )

    # Same input as every other correlation method: the abundance table as loaded.
    # SPRING's compositional model (MCLR) is applied internally by the port.
    quantitative = False
    qdat = mclr(X)
    n_lambdas = int(params.get('springNLambdas', 20) or 20)
    n_subsamples = int(params.get('springNSubsamples', 20) or 20)
    stars_threshold = float(params.get('springStarsThreshold', 0.1))
    lambda_min_ratio = float(params.get('springLambdaMinRatio', 1e-2))
    seed = int(params.get('springSeed', 17))
    subsample_ratio = float(params.get('springSubsampleFraction', 0.8))
    r_method = str(params.get('springRmethod', 'approx'))
    r_tol = float(params.get('springRtol', 1e-6))
    use_near_pd = bool(params.get('springUseNearPD', True))
    nu = float(params.get('springNu', 0.001))
    bridge_ratio = float(params.get('springRatio', 0.9))

    _progress(0.05)
    log(
        f"Running SPRING: {len(usable_samples)} samples × {len(varying_features)} taxa "
        f"(compositional MCLR, {n_lambdas} lambdas, {n_subsamples} StARS subsamples)"
    )

    def _on_progress(done: int, total: int) -> None:
        _progress(0.05 + 0.75 * (done / max(1, total)))
        log(f"  SPRING StARS progress: {done}/{total} subsamples")

    try:
        result = spring_fit(
            X,
            quantitative=quantitative,
            lambda_min_ratio=lambda_min_ratio,
            lambdaseq='data-specific',
            nlambda=n_lambdas,
            seed=seed,
            thresh=stars_threshold,
            subsample_ratio=subsample_ratio,
            rep_num=n_subsamples,
            Rtol=r_tol,
            Rmethod=r_method,
            use_nearPD=use_near_pd,
            nu=nu,
            ratio=bridge_ratio,
            verbose=False,
            validate_input=True,
        )
    except Exception as exc:
        raise RuntimeError(f'SPRING latentcor + MB/StARS failed: {exc}') from exc

    _progress(0.85)
    # The Python SpringResult exposes the selected adjacency and StARS merge
    # frequencies, but not the full latent-correlation matrix. Recompute that
    # deterministic transform with identical parameters for signed graph weights.
    try:
        latent = latentcor(
            qdat,
            types='tru',
            method=r_method,
            use_nearPD=use_near_pd,
            nu=nu,
            tol=r_tol,
            ratio=bridge_ratio,
        ).R
    except Exception as exc:
        raise RuntimeError(f'SPRING latent-correlation weight calculation failed: {exc}') from exc

    opt_index = int(result.output.stars.opt_index) - 1
    adjacency = np.asarray(result.adjacency)
    stability = np.asarray(result.output.stars.merge[opt_index], dtype=float)
    corr_matrix = {feature: {} for feature in varying_features}
    pval_matrix = {feature: {} for feature in varying_features}
    edge_count = 0
    for i, left in enumerate(varying_features):
        for j in range(i + 1, len(varying_features)):
            if int(adjacency[i, j]) == 0:
                continue
            right = varying_features[j]
            weight = float(latent[i, j])
            corr_matrix[left][right] = weight
            corr_matrix[right][left] = weight
            # The absent p-value is intentional: StARS selection is not a
            # pairwise hypothesis test.
            pval_matrix[left][right] = None
            pval_matrix[right][left] = None
            edge_count += 1

    _progress(1.0)
    diagnostics = {
        'selected_lambda': float(result.lambdaseq[opt_index]),
        'opt_index': int(result.output.stars.opt_index),
        'edge_count': edge_count,
        'n_samples': len(usable_samples),
        'n_taxa': len(varying_features),
        'mean_edge_stability': (
            float(np.mean([stability[i, j] for i in range(len(varying_features))
                           for j in range(i + 1, len(varying_features))
                           if int(adjacency[i, j]) != 0]))
            if edge_count else 0.0
        ),
        'weight_type': 'latent_correlation',
        'evidence_type': 'stars',
    }
    log(
        f"  SPRING selected lambda={diagnostics['selected_lambda']:.6g} "
        f"(index {diagnostics['opt_index']}) with {edge_count} StARS-selected edges; "
        f"mean edge stability={diagnostics['mean_edge_stability']:.4f}"
    )
    return corr_matrix, pval_matrix, diagnostics


def _compute_network_for_matrix(
    method: str,
    matrix_by_sample: dict,
    features: list,
    sample_ids: list,
    params: dict,
    progress_range: tuple[float, float] | None = None,
    drop_events: list | None = None,
    scope_label: str = '',
) -> tuple[dict, dict, dict]:
    if method == 'spieceasi':
        return _compute_spieceasi_for_matrix(
            matrix_by_sample, features, sample_ids, params, progress_range,
            drop_events=drop_events, scope_label=scope_label,
        )
    if method == 'spring':
        return _compute_spring_for_matrix(
            matrix_by_sample, features, sample_ids, params, progress_range,
            drop_events=drop_events, scope_label=scope_label,
        )
    raise ValueError(f'{method!r} is not a supported network-inference method.')


def _extract_taxon_name(feature_name: str, grouping_mode: str) -> str:
    if grouping_mode == 'all':
        return feature_name

    prefix = 'p_' if grouping_mode == 'phylum' else 'g_'
    tokens = [token.strip() for token in str(feature_name).replace(';', '|').replace(',', '|').split('|') if token.strip()]

    for token in tokens:
        lower = token.lower()
        if lower.startswith(prefix):
            cleaned = token[2:].strip('_ ').strip()
            return cleaned or f'{grouping_mode}_unknown'

    return f'{grouping_mode}_unknown'


def _group_features_before_correlation(
    sample_to_values: dict,
    ordered_features: list,
    grouping_mode: str,
    domain_label: str,
) -> tuple[dict, list]:
    if grouping_mode == 'all':
        return sample_to_values, ordered_features

    grouped_by_sample = {}
    grouped_features_seen = set()
    grouped_features = []

    for sample_id, values_by_feature in sample_to_values.items():
        grouped_values = {}
        for feature_name, value in values_by_feature.items():
            grouped_name = _extract_taxon_name(feature_name, grouping_mode)
            grouped_values[grouped_name] = grouped_values.get(grouped_name, 0.0) + float(value)
            if grouped_name not in grouped_features_seen:
                grouped_features_seen.add(grouped_name)
                grouped_features.append(grouped_name)
        grouped_by_sample[sample_id] = grouped_values

    # Ensure stable order based on original feature traversal when possible.
    ordered_grouped = []
    ordered_grouped_seen = set()
    for feature_name in ordered_features:
        grouped_name = _extract_taxon_name(feature_name, grouping_mode)
        if grouped_name not in ordered_grouped_seen:
            ordered_grouped_seen.add(grouped_name)
            ordered_grouped.append(grouped_name)
    for grouped_name in grouped_features:
        if grouped_name not in ordered_grouped_seen:
            ordered_grouped_seen.add(grouped_name)
            ordered_grouped.append(grouped_name)

    log(
        f"Applied {grouping_mode} grouping to {domain_label}: "
        f"{len(ordered_features)} original features -> {len(ordered_grouped)} grouped taxa"
    )
    return grouped_by_sample, ordered_grouped


def _restrict_features(matrix_by_sample: dict, kept_features: list) -> dict:
    """Drop every feature not in ``kept_features`` from the sample→values map.

    Grouping aggregates *all* features still present in each sample, so removing
    the filtered taxa from the matrix (not just from the feature list) is what
    actually keeps them out of the correlation step.
    """
    kept = set(kept_features)
    return {
        sample_id: {
            feature: value
            for feature, value in values.items()
            if feature in kept
        }
        for sample_id, values in matrix_by_sample.items()
    }


def _apply_taxa_filters(
    matrix_by_sample: dict,
    features: list,
    domain_label: str,
    params: dict,
    sample_totals: dict | None = None,
) -> tuple[dict, list, dict]:
    """Apply the user-configured taxa filters to one domain.

    Filters run on the individual features (OTU/ASV level) *before* taxonomic
    grouping, using each domain's own samples — the conventional order: rare and
    noisy features are removed first, then the survivors are aggregated.

    The abundance levels on which the filters operate are *relative*: every
    sample is closed to a per-sample total before prevalence/abundance/variance
    are evaluated. By default that total is this domain's own sum
    (``denominator='domain'``). In combined mode (``datasetType='both'``) the
    caller passes ``sample_totals`` with the bacteria + fungi per-sample sum, so
    mean relative abundance and variance are measured on the *same scale the
    correlation model uses* (``denominator='combined_bacteria_fungi'``). A sample
    missing from ``sample_totals`` falls back to its domain total. The input
    matrix is never transformed, because the correlation backends apply their
    own compositional transform. Returns ``(restricted_matrix, kept_features,
    report)``.

    Scope is deliberately **global**: thresholds are evaluated on every sample
    present in this domain's abundance table (all conditions), which keeps the
    candidate taxon set identical for the conditions being compared. The report
    records that scope and the denominator it used.
    """
    config = params.get('taxaFiltering') or {}
    samples = list(matrix_by_sample.keys())
    enabled = bool(config.get('enabled'))

    report = {
        'domain': domain_label,
        'enabled': enabled,
        # Thresholds see every sample in the abundance table, not only the
        # selected conditions (documented scope; see module docs above).
        'scope': 'all_samples_in_abundance_table',
        'total_samples': len(samples),
        'before': len(features),
        'after': len(features),
        'relative_abundance_denominator': (
            None if not enabled
            else ('combined_bacteria_fungi' if sample_totals is not None else 'domain')
        ),
        'removed_share_percent': 0.0,
        'aggressive_removal': False,
        'filters': {},
        'excluded': [],
    }
    if not enabled or not features or not samples:
        return matrix_by_sample, features, report

    # Per-sample totals (relative-abundance denominators). Only features of this
    # domain participate, so relative abundances are domain-relative unless the
    # caller overrides the totals (combined bacteria + fungi mode). A sample
    # whose total is zero contributes 0 to every taxon.
    domain_totals = {}
    for sample_id in samples:
        values = matrix_by_sample[sample_id]
        domain_totals[sample_id] = sum(float(values.get(feature, 0.0)) for feature in features)
    if sample_totals is None:
        resolved_totals = domain_totals
    else:
        resolved_totals = {
            sample_id: float(sample_totals.get(sample_id, default_total))
            for sample_id, default_total in domain_totals.items()
        }

    prevalence_cfg = config.get('prevalence') or {}
    mean_cfg = config.get('meanRelativeAbundance') or {}
    variability_cfg = config.get('variability') or {}

    prevalence_enabled = bool(prevalence_cfg.get('enabled'))
    mean_enabled = bool(mean_cfg.get('enabled'))
    variability_enabled = bool(variability_cfg.get('enabled'))

    prevalence_percent = float(prevalence_cfg.get('percent', 0) or 0)
    mean_threshold = float(mean_cfg.get('percent', 0) or 0) / 100.0
    variance_threshold = float(variability_cfg.get('minVariance', 0) or 0)

    removed_by = {'prevalence': 0, 'meanRelativeAbundance': 0, 'variability': 0}
    kept_features = []

    for feature in features:
        abundance_values = [float(matrix_by_sample[sample_id].get(feature, 0.0)) for sample_id in samples]
        positive_count = sum(1 for value in abundance_values if value > 0)
        prevalence_frac = positive_count / len(samples)

        relative_values = []
        for sample_id, value in zip(samples, abundance_values):
            total = resolved_totals.get(sample_id, 0.0)
            relative_values.append(value / total if total > 0 else 0.0)
        # Zero-abundance samples contribute 0 to the mean, so a taxon present in
        # few samples is diluted rather than averaged over present samples only.
        mean_relative = sum(relative_values) / len(relative_values) if relative_values else 0.0
        # Population variance of the relative abundances (metric for filter C).
        variance = (
            sum((value - mean_relative) ** 2 for value in relative_values) / len(relative_values)
            if relative_values else 0.0
        )

        reasons = []
        if prevalence_enabled and prevalence_frac < (prevalence_percent / 100.0):
            reasons.append('prevalence')
        if mean_enabled and mean_relative < mean_threshold:
            reasons.append('meanRelativeAbundance')
        if variability_enabled and variance < variance_threshold:
            reasons.append('variability')

        if reasons:
            for reason in reasons:
                removed_by[reason] += 1
            report['excluded'].append({
                'taxon': feature,
                'reasons': reasons,
                'prevalence_percent': round(prevalence_frac * 100.0, 4),
                'mean_relative_abundance_percent': round(mean_relative * 100.0, 6),
                'relative_abundance_variance': variance,
            })
        else:
            kept_features.append(feature)

    removed_total = len(features) - len(kept_features)
    report['after'] = len(kept_features)
    report['removed_share_percent'] = round((removed_total / len(features)) * 100.0, 4)
    # Guardrail: with few samples (or a combined bacteria + fungi denominator,
    # which is ~100x stricter for the minority domain) an aggressive threshold
    # can wipe out most of a domain. Surface it instead of failing later.
    report['aggressive_removal'] = bool(features) and (removed_total / len(features)) > 0.6
    report['filters'] = {
        'prevalence': {
            'enabled': prevalence_enabled,
            'threshold_percent': prevalence_percent,
            'removed': removed_by['prevalence'],
        },
        'meanRelativeAbundance': {
            'enabled': mean_enabled,
            'threshold_percent': mean_threshold * 100.0,
            'removed': removed_by['meanRelativeAbundance'],
        },
        'variability': {
            'enabled': variability_enabled,
            'threshold': variance_threshold,
            'removed': removed_by['variability'],
        },
    }

    log(
        f"Taxa filtering [{domain_label}]: {len(features)} -> {len(kept_features)} features "
        f"({len(features) - len(kept_features)} removed) across {len(samples)} samples; "
        f"relative-abundance denominator: {report['relative_abundance_denominator']}"
    )
    for name, info in report['filters'].items():
        if info['enabled']:
            log(
                f"    - {name}: removed {info['removed']} "
                f"(threshold {info.get('threshold_percent', info.get('threshold'))})"
            )
    if report['aggressive_removal']:
        log(
            f"    WARNING: the taxa filters removed {report['removed_share_percent']}% of the "
            f"{domain_label} taxa — the thresholds look very aggressive for this dataset. "
            "Consider relaxing the preset or the individual thresholds."
        )

    restricted = _restrict_features(matrix_by_sample, kept_features)
    return restricted, kept_features, report


def _build_condition_pairs(
    condition: str,
    method: str,
    params: dict,
    bacteria_by_sample: dict,
    bacteria_features: list,
    fungi_by_sample: dict,
    fungi_features: list,
    bacteria_sample_to_condition: dict,
    fungi_sample_to_condition: dict,
    progress_range: tuple[float, float] | None = None,
    drop_events: list | None = None,
) -> list:
    if method not in CORRELATION_METHODS:
        raise ValueError(
            f"Unsupported correlation method '{method}'. Use "
            f"{', '.join(sorted(CORRELATION_METHODS))}"
        )

    events_start = len(drop_events) if drop_events is not None else 0
    # Cross-domain pairs are fitted on one combined bac+fun matrix.
    scope_label = f'bac+fun · {condition}'

    bacteria_condition_samples = {
        sample_id
        for sample_id, sample_condition in bacteria_sample_to_condition.items()
        if sample_condition == condition
    }
    fungi_condition_samples = {
        sample_id
        for sample_id, sample_condition in fungi_sample_to_condition.items()
        if sample_condition == condition
    }

    common_samples = sorted(
        bacteria_condition_samples
        & fungi_condition_samples
        & set(bacteria_by_sample.keys())
        & set(fungi_by_sample.keys())
    )

    log(
        f"  Condition '{condition}' cross-domain matrix: {len(common_samples)} shared samples × "
        f"{len(bacteria_features)} bacteria features × {len(fungi_features)} fungi features"
    )

    rows = []
    corr_matrix = {}
    pval_matrix = {}
    diagnostics = {}
    if method in MATRIX_METHODS:
        combined_features = [f'bac::{name}' for name in bacteria_features] + [f'fun::{name}' for name in fungi_features]
        combined_matrix_by_sample = {}
        for sample_id in common_samples:
            sample_values = {}
            for feature in bacteria_features:
                sample_values[f'bac::{feature}'] = bacteria_by_sample[sample_id].get(feature, 0.0)
            for feature in fungi_features:
                sample_values[f'fun::{feature}'] = fungi_by_sample[sample_id].get(feature, 0.0)
            combined_matrix_by_sample[sample_id] = sample_values
        if method == 'sparcc':
            corr_matrix, pval_matrix = _compute_sparcc_for_matrix(
                combined_matrix_by_sample,
                combined_features,
                common_samples,
                params,
                progress_range,
                drop_events=drop_events,
                scope_label=scope_label,
            )
        elif method == 'cclasso':
            corr_matrix, pval_matrix, diagnostics = _compute_cclasso_for_matrix(
                combined_matrix_by_sample,
                combined_features,
                common_samples,
                params,
                progress_range,
                drop_events=drop_events,
                scope_label=scope_label,
            )
        else:
            corr_matrix, pval_matrix, diagnostics = _compute_network_for_matrix(
                method,
                combined_matrix_by_sample,
                combined_features,
                common_samples,
                params,
                progress_range,
                drop_events=drop_events,
                scope_label=scope_label,
            )
            if diagnostics.get('edge_count', 0) > 0:
                log(
                    f"  {method.upper()} cross-domain selected lambda="
                    f"{diagnostics['selected_lambda']:.6g} with {diagnostics['edge_count']} edges"
                )
        if drop_events is not None:
            for event in drop_events[events_start:]:
                event.setdefault('condition', condition)

    for bac in bacteria_features:
        for fun in fungi_features:
            x_values = []
            y_values = []
            for sample_id in common_samples:
                bac_value = bacteria_by_sample[sample_id].get(bac)
                fun_value = fungi_by_sample[sample_id].get(fun)
                if bac_value is None or fun_value is None:
                    continue
                x_values.append(bac_value)
                y_values.append(fun_value)

            if len(x_values) < 3:
                continue

            if method in MATRIX_METHODS:
                correlation, p_value = _get_pair_stat(
                    corr_matrix,
                    pval_matrix,
                    f'bac::{bac}',
                    f'fun::{fun}',
                )
            else:
                correlation, p_value = _compute_spearman(x_values, y_values)
            if math.isnan(correlation) or (
                p_value is not None and not math.isfinite(p_value)
            ):
                continue

            row = {
                'row': bac,
                'column': fun,
                'cor': correlation,
                'p': p_value,
                'method': method,
                'condition': condition,
                'row_kind': 'bacteria',
                'column_kind': 'fungi',
                'weight_type': diagnostics.get('weight_type', 'correlation'),
                'evidence_type': diagnostics.get('evidence_type', 'p_value'),
                'method_diagnostics': diagnostics,
            }
            if method in NETWORK_METHODS:
                row['selected'] = True
                row['selection_method'] = 'StARS'
            rows.append(row)

    return rows


def _build_both_domain_pairs(
    condition: str,
    method: str,
    params: dict,
    bacteria_by_sample: dict,
    bacteria_features: list,
    fungi_by_sample: dict,
    fungi_features: list,
    bacteria_sample_to_condition: dict,
    fungi_sample_to_condition: dict,
    progress_range: tuple[float, float] | None = None,
    drop_events: list | None = None,
) -> list:
    events_start = len(drop_events) if drop_events is not None else 0
    if method in NETWORK_METHODS:
        bacteria_condition_samples = {
            sample_id
            for sample_id, sample_condition in bacteria_sample_to_condition.items()
            if sample_condition == condition
        }
        fungi_condition_samples = {
            sample_id
            for sample_id, sample_condition in fungi_sample_to_condition.items()
            if sample_condition == condition
        }
        common_samples = sorted(
            bacteria_condition_samples
            & fungi_condition_samples
            & set(bacteria_by_sample.keys())
            & set(fungi_by_sample.keys())
        )

        combined_features = [f'bac::{name}' for name in bacteria_features] + [f'fun::{name}' for name in fungi_features]
        combined_matrix_by_sample = {}
        for sample_id in common_samples:
            sample_values = {}
            for feature in bacteria_features:
                sample_values[f'bac::{feature}'] = bacteria_by_sample[sample_id].get(feature, 0.0)
            for feature in fungi_features:
                sample_values[f'fun::{feature}'] = fungi_by_sample[sample_id].get(feature, 0.0)
            combined_matrix_by_sample[sample_id] = sample_values

        corr_matrix, pval_matrix, diagnostics = _compute_network_for_matrix(
            method,
            combined_matrix_by_sample,
            combined_features,
            common_samples,
            params,
            progress_range,
            drop_events=drop_events,
            scope_label=f'bac+fun · {condition}',
        )

        rows = []
        for i, left in enumerate(combined_features):
            for j in range(i + 1, len(combined_features)):
                right = combined_features[j]
                correlation, p_value = _get_pair_stat(corr_matrix, pval_matrix, left, right)
                if math.isnan(correlation) or (
                    p_value is not None and not math.isfinite(p_value)
                ):
                    continue
                left_kind, left_name = ('bacteria', left[5:]) if left.startswith('bac::') else ('fungi', left[5:])
                right_kind, right_name = ('bacteria', right[5:]) if right.startswith('bac::') else ('fungi', right[5:])
                row = {
                    'row': left_name,
                    'column': right_name,
                    'cor': correlation,
                    'p': p_value,
                    'method': method,
                    'condition': condition,
                    'row_kind': left_kind,
                    'column_kind': right_kind,
                    'weight_type': diagnostics.get('weight_type', 'latent_correlation'),
                    'evidence_type': 'stars',
                    'selected': True,
                    'selection_method': 'StARS',
                    'method_diagnostics': diagnostics,
                }
                rows.append(row)

        log(
            f"  Combined pairs ({method.upper()}): total={len(rows)}, selected_edges="
            f"{diagnostics.get('edge_count', 0)}, selected_lambda={diagnostics.get('selected_lambda', math.nan):.6g}"
        )
        if drop_events is not None:
            for event in drop_events[events_start:]:
                event.setdefault('condition', condition)
        return rows

    # In combined mode, include all biologically relevant edges:
    # 1) bacteria-bacteria, 2) fungi-fungi, 3) bacteria-fungi.
    # Split the progress window evenly across the three sub-matrices.
    lo, hi = progress_range or (0.0, 0.0)
    third = (hi - lo) / 3.0
    bac_pairs = _build_same_domain_pairs(
        condition,
        method,
        params,
        bacteria_by_sample,
        bacteria_features,
        bacteria_sample_to_condition,
        'bacteria',
        progress_range=(lo, lo + third),
        drop_events=drop_events,
    )
    fun_pairs = _build_same_domain_pairs(
        condition,
        method,
        params,
        fungi_by_sample,
        fungi_features,
        fungi_sample_to_condition,
        'fungi',
        progress_range=(lo + third, lo + 2 * third),
        drop_events=drop_events,
    )
    cross_pairs = _build_condition_pairs(
        condition,
        method,
        params,
        bacteria_by_sample,
        bacteria_features,
        fungi_by_sample,
        fungi_features,
        bacteria_sample_to_condition,
        fungi_sample_to_condition,
        progress_range=(lo + 2 * third, hi),
        drop_events=drop_events,
    )
    log(
        f"  Combined pairs: bacteria-bacteria={len(bac_pairs)}, "
        f"fungi-fungi={len(fun_pairs)}, bacteria-fungi={len(cross_pairs)}"
    )
    return bac_pairs + fun_pairs + cross_pairs


def _build_same_domain_pairs(
    condition: str,
    method: str,
    params: dict,
    domain_by_sample: dict,
    domain_features: list,
    sample_to_condition: dict,
    kind_label: str,
    progress_range: tuple[float, float] | None = None,
    drop_events: list | None = None,
) -> list:
    if method not in CORRELATION_METHODS:
        raise ValueError(
            f"Unsupported correlation method '{method}'. Use "
            f"{', '.join(sorted(CORRELATION_METHODS))}"
        )

    events_start = len(drop_events) if drop_events is not None else 0
    scope_label = f'{kind_label} · {condition}'

    condition_samples = {
        sample_id
        for sample_id, sample_condition in sample_to_condition.items()
        if sample_condition == condition
    }

    common_samples = sorted(condition_samples & set(domain_by_sample.keys()))
    log(
        f"  Condition '{condition}' {kind_label} matrix: {len(common_samples)} samples × "
        f"{len(domain_features)} features → {len(common_samples) * len(domain_features)} data points"
    )
    
    max_pairs = len(domain_features) * (len(domain_features) - 1) // 2
    log(f"  Max possible pairs: {max_pairs} (combinations of {len(domain_features)} features)")
    rows = []
    corr_matrix = {}
    pval_matrix = {}
    diagnostics = {}
    if method in MATRIX_METHODS:
        if method == 'sparcc':
            corr_matrix, pval_matrix = _compute_sparcc_for_matrix(
                domain_by_sample,
                domain_features,
                common_samples,
                params,
                progress_range,
                drop_events=drop_events,
                scope_label=scope_label,
            )
        elif method == 'cclasso':
            corr_matrix, pval_matrix, diagnostics = _compute_cclasso_for_matrix(
                domain_by_sample,
                domain_features,
                common_samples,
                params,
                progress_range,
                drop_events=drop_events,
                scope_label=scope_label,
            )
        else:
            corr_matrix, pval_matrix, diagnostics = _compute_network_for_matrix(
                method,
                domain_by_sample,
                domain_features,
                common_samples,
                params,
                progress_range,
                drop_events=drop_events,
                scope_label=scope_label,
            )
            if diagnostics.get('edge_count', 0) > 0:
                log(
                    f"  {method.upper()} {kind_label} selected lambda="
                    f"{diagnostics['selected_lambda']:.6g} with {diagnostics['edge_count']} edges"
                )
        if drop_events is not None:
            for event in drop_events[events_start:]:
                event.setdefault('condition', condition)

    for i in range(len(domain_features)):
        left = domain_features[i]
        for j in range(i + 1, len(domain_features)):
            right = domain_features[j]
            x_values = []
            y_values = []
            for sample_id in common_samples:
                left_value = domain_by_sample[sample_id].get(left)
                right_value = domain_by_sample[sample_id].get(right)
                if left_value is None or right_value is None:
                    continue
                x_values.append(left_value)
                y_values.append(right_value)

            if len(x_values) < 3:
                continue

            if method in MATRIX_METHODS:
                correlation, p_value = _get_pair_stat(corr_matrix, pval_matrix, left, right)
            else:
                correlation, p_value = _compute_spearman(x_values, y_values)
            if math.isnan(correlation) or (
                p_value is not None and not math.isfinite(p_value)
            ):
                continue

            row = {
                'row': left,
                'column': right,
                'cor': correlation,
                'p': p_value,
                'method': method,
                'condition': condition,
                'row_kind': kind_label,
                'column_kind': kind_label,
                'weight_type': diagnostics.get('weight_type', 'correlation'),
                'evidence_type': diagnostics.get('evidence_type', 'p_value'),
                'method_diagnostics': diagnostics,
            }
            if method in NETWORK_METHODS:
                row['selected'] = True
                row['selection_method'] = 'StARS'
            rows.append(row)

    log(
        f"  Pairs computed: {len(rows)} ({method} correlation on {len(common_samples)} samples per pair, "
        f"minimum 3 samples required)"
    )
    return rows


def _write_pairs_csv(job_id: str, condition: str, rows: list) -> str:
    file_name = f"corr_pairs_{job_id}_{_slug(condition)}.csv"
    output_path = os.path.join(tempfile.gettempdir(), file_name)
    fieldnames = [
        'row', 'column', 'cor', 'p', 'weight_type', 'evidence_type',
        'selection_method', 'selected',
    ]
    with open(output_path, 'w', newline='', encoding='utf-8') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                'row': row['row'],
                'column': row['column'],
                'cor': row['cor'],
                'p': row['p'],
                'weight_type': row.get('weight_type', 'correlation'),
                'evidence_type': row.get('evidence_type', 'p_value'),
                'selection_method': row.get('selection_method', ''),
                'selected': row.get('selected', ''),
            })
    return output_path


def _build_graph_from_pairs(
    rows: list,
    p_value_threshold: float | None,
    overabundance_by_node: dict | None = None,
) -> dict:
    if p_value_threshold is None:
        selected_rows = [row for row in rows if row.get('selected', False)]
    else:
        selected_rows = [
            row for row in rows
            if row.get('p') is not None and row['p'] <= p_value_threshold
        ]
    node_ids = set()
    node_kinds = {}
    edges = []
    for row in selected_rows:
        source = row['row']
        target = row['column']
        node_ids.add(source)
        node_ids.add(target)
        node_kinds[source] = row.get('row_kind', 'feature')
        node_kinds[target] = row.get('column_kind', 'feature')
        edges.append({
            'source': source,
            'target': target,
            'weight': row['cor'],
            'p_value': row['p'],
            'method': row['method'],
            'weight_type': row.get('weight_type', 'correlation'),
            'evidence_type': row.get('evidence_type', 'p_value'),
            'selection_method': row.get('selection_method'),
        })

    nodes = []
    overabundance_by_node = overabundance_by_node or {}
    for node_id in sorted(node_ids):
        kind = node_kinds.get(node_id, 'feature')
        nodes.append({
            'id': node_id,
            'label': node_id,
            'kind': kind,
            'overabundance': float(overabundance_by_node.get((node_id, kind), 1.0)),
        })

    stats = {
        'total_pairs': len(rows),
        'node_count': len(nodes),
        'edge_count': len(edges),
        'p_value_threshold': p_value_threshold,
        'p_value_threshold_applied': p_value_threshold is not None,
        'weight_type': rows[0].get('weight_type', 'correlation') if rows else 'correlation',
        'evidence_type': rows[0].get('evidence_type', 'p_value') if rows else 'p_value',
        'method_diagnostics': rows[0].get('method_diagnostics', {}) if rows else {},
    }
    if p_value_threshold is None:
        stats['selected_pairs'] = len(selected_rows)
    else:
        stats['significant_pairs'] = len(selected_rows)

    return {
        'nodes': nodes,
        'edges': edges,
        'stats': stats,
    }


def _condition_sample_stats(
    condition: str,
    sample_to_condition: dict,
    bacteria_by_sample: dict,
    fungi_by_sample: dict,
    dataset_type: str,
) -> dict:
    condition_samples = {
        sample_id
        for sample_id, sample_condition in sample_to_condition.items()
        if sample_condition == condition
    }

    stats = {
        'condition_samples_in_metadata': len(condition_samples),
        'condition_samples_after_dataset_filter': 0,
    }

    if dataset_type == 'bacteria':
        stats['condition_samples_after_dataset_filter'] = len(condition_samples & set(bacteria_by_sample.keys()))
    elif dataset_type == 'fungi':
        stats['condition_samples_after_dataset_filter'] = len(condition_samples & set(fungi_by_sample.keys()))
    else:
        stats['condition_samples_after_dataset_filter'] = len(
            condition_samples
            & set(bacteria_by_sample.keys())
            & set(fungi_by_sample.keys())
        )

    return stats


def _mean_by_feature_for_samples(sample_ids: set[str], by_sample: dict, features: list) -> dict:
    means = {}
    for feature in features:
        values = []
        for sample_id in sample_ids:
            sample_values = by_sample.get(sample_id)
            if not sample_values:
                continue
            value = sample_values.get(feature)
            if value is None:
                continue
            values.append(float(value))
        means[feature] = (sum(values) / len(values)) if values else 0.0
    return means


def _compute_overabundance_map(
    condition: str,
    dataset_type: str,
    bacteria_by_sample: dict,
    bacteria_features: list,
    bacteria_sample_to_condition: dict,
    fungi_by_sample: dict,
    fungi_features: list,
    fungi_sample_to_condition: dict,
) -> dict:
    overabundance = {}

    if dataset_type in {'bacteria', 'both'}:
        condition_samples = {
            sample_id
            for sample_id, sample_condition in bacteria_sample_to_condition.items()
            if sample_condition == condition and sample_id in bacteria_by_sample
        }
        overall_samples = set(bacteria_by_sample.keys())
        condition_means = _mean_by_feature_for_samples(condition_samples, bacteria_by_sample, bacteria_features)
        overall_means = _mean_by_feature_for_samples(overall_samples, bacteria_by_sample, bacteria_features)

        for feature in bacteria_features:
            overall = overall_means.get(feature, 0.0)
            condition_mean = condition_means.get(feature, 0.0)
            ratio = condition_mean / overall if overall > 0 else 1.0
            overabundance[(feature, 'bacteria')] = float(ratio)

    if dataset_type in {'fungi', 'both'}:
        condition_samples = {
            sample_id
            for sample_id, sample_condition in fungi_sample_to_condition.items()
            if sample_condition == condition and sample_id in fungi_by_sample
        }
        overall_samples = set(fungi_by_sample.keys())
        condition_means = _mean_by_feature_for_samples(condition_samples, fungi_by_sample, fungi_features)
        overall_means = _mean_by_feature_for_samples(overall_samples, fungi_by_sample, fungi_features)

        for feature in fungi_features:
            overall = overall_means.get(feature, 0.0)
            condition_mean = condition_means.get(feature, 0.0)
            ratio = condition_mean / overall if overall > 0 else 1.0
            overabundance[(feature, 'fungi')] = float(ratio)

    return overabundance


def _condition_sample_stats_both(
    condition: str,
    bacteria_sample_to_condition: dict,
    fungi_sample_to_condition: dict,
    bacteria_by_sample: dict,
    fungi_by_sample: dict,
) -> dict:
    bacteria_condition_samples = {
        sample_id
        for sample_id, sample_condition in bacteria_sample_to_condition.items()
        if sample_condition == condition
    }
    fungi_condition_samples = {
        sample_id
        for sample_id, sample_condition in fungi_sample_to_condition.items()
        if sample_condition == condition
    }

    bacteria_after_filter = bacteria_condition_samples & set(bacteria_by_sample.keys())
    fungi_after_filter = fungi_condition_samples & set(fungi_by_sample.keys())
    cross_after_filter = bacteria_after_filter & fungi_after_filter

    return {
        'condition_samples_in_metadata': len(bacteria_condition_samples | fungi_condition_samples),
        'condition_samples_after_dataset_filter': len(bacteria_after_filter | fungi_after_filter),
        'condition_bacteria_samples_after_filter': len(bacteria_after_filter),
        'condition_fungi_samples_after_filter': len(fungi_after_filter),
        'condition_cross_domain_samples_after_filter': len(cross_after_filter),
    }

def run(job_id: str, files: dict, params: dict) -> dict:
    """Core analysis logic — replace with your own."""
    bacteria_files = files.get("bacteria", [])
    fungi_files = files.get("fungi", [])
    bacteria_metadata_files = files.get("bacteria_metadata", [])
    fungi_metadata_files = files.get("fungi_metadata", [])

    dataset_type = params.get("datasetType")
    correlation_method = str(params.get("correlationMethod") or 'spearman').lower()
    condition_column = str(params.get("conditionColumn") or 'Combination_treat2').strip() or 'Combination_treat2'
    condition_a = params.get("conditionA")
    condition_b = params.get("conditionB")
    grouping_mode = params.get('groupingMode', 'all')
    p_value_threshold = (
        float(params.get('pValueThreshold', 0.05))
        if correlation_method in P_VALUE_METHODS else None
    )
    # Minimum absolute edge score. It is a correlation threshold for pairwise
    # methods, StARS stability for SPIEC-EASI, and latent correlation for SPRING.
    correlation_threshold = abs(float(params.get('correlationThreshold', 0) or 0))

    log(f"Starting analysis for job {job_id}")
    log(f"Dataset type: {dataset_type}")
    log(f"Correlation method: {correlation_method}")
    log(f"Grouping mode: {grouping_mode}")
    if correlation_method == 'spieceasi':
        log(f"StARS edge-stability filter (minimum): {correlation_threshold}")
    elif correlation_method == 'spring':
        log(f"Correlation filter (minimum |r|, latent correlation): {correlation_threshold}")
    else:
        log(f"Correlation filter (minimum |r|): {correlation_threshold}")
    log(f"Condition column: {condition_column}")
    log(f"Condition comparison: {condition_a} vs {condition_b}")
    log(f"Bacteria files: {bacteria_files}")
    log(f"Fungi files: {fungi_files}")
    log(f"Bacteria metadata files: {bacteria_metadata_files}")
    log(f"Fungi metadata files: {fungi_metadata_files}")

    if not bacteria_metadata_files and not fungi_metadata_files:
        raise ValueError('At least one metadata file (bacteria or fungi) is required')
    if dataset_type not in {'bacteria', 'fungi', 'both'}:
        raise ValueError("datasetType must be one of: bacteria, fungi, both")

    if dataset_type in {'bacteria', 'both'} and not bacteria_files:
        raise ValueError('Bacteria files are required for the selected dataset type')
    if dataset_type in {'fungi', 'both'} and not fungi_files:
        raise ValueError('Fungi files are required for the selected dataset type')
    if not condition_a or not condition_b:
        raise ValueError('Both conditionA and conditionB are required')
    if condition_a == condition_b:
        raise ValueError('conditionA and conditionB must be different')
    if grouping_mode not in {'all', 'phylum', 'genus'}:
        raise ValueError("groupingMode must be one of: all, phylum, genus")
    if correlation_method not in CORRELATION_METHODS:
        raise ValueError(
            f"correlationMethod must be one of: {', '.join(sorted(CORRELATION_METHODS))}"
        )
    if correlation_threshold < 0 or correlation_threshold > 1:
        raise ValueError("correlationThreshold must be between 0 and 1")

    # Load metadata separately for each domain
    bacteria_sample_to_condition = _load_metadata_conditions(
        bacteria_metadata_files,
        'Bacteria',
        condition_column,
    ) if bacteria_metadata_files else {}
    fungi_sample_to_condition = _load_metadata_conditions(
        fungi_metadata_files,
        'Fungi',
        condition_column,
    ) if fungi_metadata_files else {}

    # Combine metadata for loading abundance matrices (for transpose detection)
    all_known_sample_ids = set(bacteria_sample_to_condition.keys()) | set(fungi_sample_to_condition.keys())
    bacteria_by_sample, bacteria_features = _load_abundance_matrices(
        bacteria_files,
        all_known_sample_ids
    ) if bacteria_files else ({}, [])
    fungi_by_sample, fungi_features = _load_abundance_matrices(
        fungi_files,
        all_known_sample_ids
    ) if fungi_files else ({}, [])

    # ── Taxa filtering (per domain, before grouping) ──────────────────────────
    # Applied here so every downstream step (grouping, correlation, graphs) sees
    # only the surviving taxa. Disabled master toggle ⇒ both calls are no-ops.
    # Scope is deliberately global: thresholds use every sample of the domain's
    # abundance table, so both selected conditions start from the same candidate
    # taxon set (and the same graph node pool).
    taxa_filtering_config = params.get('taxaFiltering') or {}
    log(f"Taxa filtering enabled: {bool(taxa_filtering_config.get('enabled'))}")
    log(
        "Taxa filtering scope: every sample of the abundance table (all conditions) — "
        "both selected conditions start from the same candidate taxon set"
    )
    taxa_filtering_report = {}

    # Relative-abundance denominator for filters B (mean) and C (variance). In
    # combined mode the correlation methods close the fused bac::/fun:: matrix on
    # bacteria + fungi *together*, so the filters must measure abundance on that
    # same scale: with a per-domain denominator the minority domain looks ~100x
    # more abundant to the filter than it is to the model (variance ~10^4x).
    # Prevalence is scale-free and is unaffected by this choice.
    filter_totals: dict | None = None
    if dataset_type == 'both' and bacteria_by_sample and fungi_by_sample:
        bac_totals = {
            sample_id: sum(float(values.get(f, 0.0)) for f in bacteria_features)
            for sample_id, values in bacteria_by_sample.items()
        }
        fun_totals = {
            sample_id: sum(float(values.get(f, 0.0)) for f in fungi_features)
            for sample_id, values in fungi_by_sample.items()
        }
        filter_totals = {
            sample_id: bac_totals.get(sample_id, 0.0) + fun_totals.get(sample_id, 0.0)
            for sample_id in set(bac_totals) | set(fun_totals)
        }
        log(
            "Taxa-filtering denominator: combined bacteria + fungi per-sample total "
            "(same scale as the combined matrix fed to the correlation model)"
        )

    taxa_filtering_report['meta'] = {
        'scope': 'all_samples_in_abundance_table',
        'relative_abundance_denominator': (
            'combined_bacteria_fungi' if filter_totals is not None else 'domain'
        ),
        'method_removals': [],
        'method_removals_note': (
            'Taxa/samples discarded by the correlation method itself (zero variance in the '
            'analysed sample subset, all-zero samples) after the user filters above. Scopes '
            'are per condition and matrix, so counts may overlap when datasetType="both"; the '
            'exact model input size is graph.stats.method_diagnostics.n_taxa.'
        ),
    }
    if bacteria_by_sample:
        bacteria_by_sample, bacteria_features, bac_filter_report = _apply_taxa_filters(
            bacteria_by_sample,
            bacteria_features,
            'bacteria',
            params,
            sample_totals=filter_totals,
        )
        taxa_filtering_report['bacteria'] = bac_filter_report
    if fungi_by_sample:
        fungi_by_sample, fungi_features, fun_filter_report = _apply_taxa_filters(
            fungi_by_sample,
            fungi_features,
            'fungi',
            params,
            sample_totals=filter_totals,
        )
        taxa_filtering_report['fungi'] = fun_filter_report

    # A domain that feeds the correlation must keep enough features to be
    # meaningful; otherwise fail loudly instead of producing an empty graph.
    if dataset_type in {'bacteria', 'both'} and bacteria_by_sample and len(bacteria_features) < 2:
        raise ValueError(
            f"Taxa filtering left only {len(bacteria_features)} bacteria feature(s). "
            "Disable or relax the filters to keep at least 2 features."
        )
    if dataset_type in {'fungi', 'both'} and fungi_by_sample and len(fungi_features) < 2:
        raise ValueError(
            f"Taxa filtering left only {len(fungi_features)} fungi feature(s). "
            "Disable or relax the filters to keep at least 2 features."
        )

    bacteria_by_sample, bacteria_features = _group_features_before_correlation(
        bacteria_by_sample,
        bacteria_features,
        grouping_mode,
        'bacteria',
    ) if bacteria_by_sample else ({}, [])
    fungi_by_sample, fungi_features = _group_features_before_correlation(
        fungi_by_sample,
        fungi_features,
        grouping_mode,
        'fungi',
    ) if fungi_by_sample else ({}, [])
    
    log("Loaded abundance matrices:")
    log(f"  Bacteria: {len(bacteria_by_sample)} samples × {len(bacteria_features)} features")
    log(f"  Fungi: {len(fungi_by_sample)} samples × {len(fungi_features)} features")

    if bacteria_by_sample:
        bacteria_overlap = len(set(bacteria_by_sample.keys()) & set(bacteria_sample_to_condition.keys()))
        bacteria_missing = len(set(bacteria_sample_to_condition.keys()) - set(bacteria_by_sample.keys()))
        log(f"  Bacteria/bacteria_metadata sample overlap: {bacteria_overlap}/{len(bacteria_sample_to_condition)} "
            f"(missing: {bacteria_missing})")
        if bacteria_overlap == 0 and dataset_type in {'bacteria', 'both'}:
            log("  WARNING: No overlap between bacteria and bacteria metadata samples!")
    if fungi_by_sample:
        fungi_overlap = len(set(fungi_by_sample.keys()) & set(fungi_sample_to_condition.keys()))
        fungi_missing = len(set(fungi_sample_to_condition.keys()) - set(fungi_by_sample.keys()))
        log(f"  Fungi/fungi_metadata sample overlap: {fungi_overlap}/{len(fungi_sample_to_condition)} "
            f"(missing: {fungi_missing})")
        if fungi_overlap == 0 and dataset_type in {'fungi', 'both'}:
            log("  WARNING: No overlap between fungi and fungi metadata samples!")

    # Build graph artifacts for selected conditions only.
    selected_conditions = [condition_a, condition_b]
    condition_graphs = []
    # Method-side removals (constant taxa / empty samples) collected per run and
    # attached to the taxa-filtering report below — see _record_method_drop.
    method_drop_events: list = []

    for idx, condition in enumerate(selected_conditions):
        if not condition:
            continue
        report_progress(10 + idx * 30)
        
        # Use domain-specific metadata for filtering.
        if dataset_type == 'bacteria':
            sample_to_condition = bacteria_sample_to_condition
            sample_stats = _condition_sample_stats(
                condition,
                sample_to_condition,
                bacteria_by_sample,
                fungi_by_sample,
                dataset_type,
            )
        elif dataset_type == 'fungi':
            sample_to_condition = fungi_sample_to_condition
            sample_stats = _condition_sample_stats(
                condition,
                sample_to_condition,
                bacteria_by_sample,
                fungi_by_sample,
                dataset_type,
            )
        else:
            sample_stats = _condition_sample_stats_both(
                condition,
                bacteria_sample_to_condition,
                fungi_sample_to_condition,
                bacteria_by_sample,
                fungi_by_sample,
            )
            sample_to_condition = {**fungi_sample_to_condition, **bacteria_sample_to_condition}
        log(
            f"Condition '{condition}': {sample_stats['condition_samples_in_metadata']} samples in metadata, "
            f"{sample_stats['condition_samples_after_dataset_filter']} after dataset filtering"
        )
        log(f"Building correlation pair table for condition '{condition}'")
        method = correlation_method
        
        # Calculate condition samples for this condition
        condition_samples = {
            sample_id
            for sample_id, sample_condition in sample_to_condition.items()
            if sample_condition == condition
        }
        
        if dataset_type == 'both':
            bac_in_cond = len(set(bacteria_by_sample.keys()) & condition_samples)
            fun_in_cond = len(set(fungi_by_sample.keys()) & condition_samples)
            log(f"  Available bacteria samples: {bac_in_cond} | Available fungi samples: {fun_in_cond}")
            all_pairs = _build_both_domain_pairs(
                condition,
                method,
                params,
                bacteria_by_sample,
                bacteria_features,
                fungi_by_sample,
                fungi_features,
                bacteria_sample_to_condition,
                fungi_sample_to_condition,
                progress_range=(10 + idx * 30, 25 + idx * 30),
                drop_events=method_drop_events,
            )
        elif dataset_type == 'bacteria':
            all_pairs = _build_same_domain_pairs(
                condition,
                method,
                params,
                bacteria_by_sample,
                bacteria_features,
                sample_to_condition,
                'bacteria',
                progress_range=(10 + idx * 30, 25 + idx * 30),
                drop_events=method_drop_events,
            )
        else:
            all_pairs = _build_same_domain_pairs(
                condition,
                method,
                params,
                fungi_by_sample,
                fungi_features,
                sample_to_condition,
                'fungi',
                progress_range=(10 + idx * 30, 25 + idx * 30),
                drop_events=method_drop_events,
            )

        report_progress(25 + idx * 30)
        if p_value_threshold is None:
            log(
                f"Filtering StARS-selected network edges for condition '{condition}' "
                f"(minimum score={correlation_threshold}; no pairwise p-values are produced)"
            )
            significant_pairs = [row for row in all_pairs if row.get('selected', False)]
        else:
            log(
                f"Filtering pairs for condition '{condition}' "
                f"(p <= {p_value_threshold}, |r| >= {correlation_threshold})"
            )
            significant_pairs = [
                row for row in all_pairs
                if row.get('p') is not None and row['p'] <= p_value_threshold
            ]
        filtered_pairs = [
            row for row in significant_pairs
            if abs(float(row['cor'])) >= correlation_threshold
        ]
        log(
            f"  All pairs computed: {len(all_pairs)} | "
            f"{'Significant' if p_value_threshold is not None else 'Selected'} "
            f"({'p ≤ ' + str(p_value_threshold) if p_value_threshold is not None else 'StARS'}): "
            f"{len(significant_pairs)} | "
            f"After |r| ≥ {correlation_threshold}: {len(filtered_pairs)}"
        )
        pair_file_path = _write_pairs_csv(job_id, condition, filtered_pairs)
        overabundance_map = _compute_overabundance_map(
            condition,
            dataset_type,
            bacteria_by_sample,
            bacteria_features,
            bacteria_sample_to_condition,
            fungi_by_sample,
            fungi_features,
            fungi_sample_to_condition,
        )

        graph = _build_graph_from_pairs(filtered_pairs, p_value_threshold, overabundance_map)
        graph['stats']['correlation_threshold'] = correlation_threshold
        graph['stats']['weight_type'] = (
            'selection_stability' if correlation_method == 'spieceasi'
            else 'latent_correlation' if correlation_method == 'spring'
            else 'correlation'
        )
        graph['stats']['evidence_type'] = (
            'stars' if correlation_method in NETWORK_METHODS
            else 'bootstrap_p_value' if correlation_method == 'cclasso'
            else 'permutation_p_value' if correlation_method == 'sparcc'
            else 'pairwise_p_value'
        )
        graph['stats']['correlation_threshold_label'] = (
            'minimum_stars_edge_stability'
            if correlation_method == 'spieceasi'
            else 'minimum_absolute_correlation'
        )
        graph_id = f"{job_id}:{_slug(condition)}"
        graph['stats'].update(sample_stats)
        condition_graphs.append({
            'graph_id': graph_id,
            'condition': condition,
            'correlation_method': correlation_method,
            'p_value_threshold_applied': p_value_threshold is not None,
            'evidence_type': (
                'stars'
                if correlation_method in NETWORK_METHODS
                else 'bootstrap_p_value'
                if correlation_method == 'cclasso'
                else 'permutation_p_value'
                if correlation_method == 'sparcc'
                else 'pairwise_p_value'
            ),
            'weight_type': (
                'selection_stability'
                if correlation_method == 'spieceasi'
                else 'latent_correlation'
                if correlation_method == 'spring'
                else 'correlation'
            ),
            'pair_file': {
                'format': 'csv',
                'path': pair_file_path,
                'row_count': len(filtered_pairs),
            },
            'graph': graph,
            'clustering': {
                'status': 'not_requested',
                'algorithm': None,
                'seed_node': None,
                'clusters': [],
                'message': 'Clustering algorithm not implemented yet.',
            },
        })

    # ── Method-side removals → taxa-filtering report ─────────────────────────
    # These used to exist only as log lines, which made report before/after
    # disagree with the number of taxa the model actually saw.
    taxa_filtering_report['meta']['method_removals'] = method_drop_events
    for event in method_drop_events:
        log(
            f"Method-side removal [{event['scope']} · {event['method']}]: "
            f"{event['count']} {event['item_kind']} ({event['kind']})"
        )

    # Simulate remaining pipeline work and complete progress.
    for step, pct in enumerate([85, 100]):
        time.sleep(0.5)  # replace with real work
        report_progress(pct)
        log(f"Finalization step {step + 1} complete")

    # ── Build result ──────────────────────────────────
    result = {
        "summary": {
            "n_bacteria_files": len(bacteria_files),
            "n_fungi_files": len(fungi_files),
            "n_bacteria_metadata_files": len(bacteria_metadata_files),
            "n_fungi_metadata_files": len(fungi_metadata_files),
            "dataset_type": dataset_type,
            "correlation_method": correlation_method,
            "method_settings": json.dumps(
                _method_settings(correlation_method, params),
                sort_keys=True,
                separators=(',', ':'),
            ),
            "condition_a": condition_a,
            "condition_b": condition_b,
            "selected_conditions": [c for c in selected_conditions if c],
            "grouping_mode": grouping_mode,
            "p_value_threshold": p_value_threshold,
            "p_value_threshold_applied": p_value_threshold is not None,
            "evidence_type": (
                'stars' if correlation_method in NETWORK_METHODS
                else 'bootstrap_p_value' if correlation_method == 'cclasso'
                else 'permutation_p_value' if correlation_method == 'sparcc'
                else 'pairwise_p_value'
            ),
            "weight_type": (
                'selection_stability' if correlation_method == 'spieceasi'
                else 'latent_correlation' if correlation_method == 'spring'
                else 'correlation'
            ),
            "correlation_threshold": correlation_threshold,
            "taxa_filtering_enabled": bool(taxa_filtering_config.get('enabled')),
            "taxa_filtering": taxa_filtering_report,
        },
        # These are reserved for future, real analysis outputs. Never populate
        # them with illustrative values in a statistical-results payload.
        "charts": [],
        "table": [],
        "graph_pipeline": {
            "version": 1,
            "status": "completed",
            "condition_graphs": condition_graphs,
        },
    }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", required=True, help="Path to params JSON file")
    parser.add_argument(
        "--output",
        default=None,
        help="Where to write the result JSON. Defaults to the OS temporary directory.",
    )
    args = parser.parse_args()

    with open(args.params, "r", encoding="utf-8") as f:
        payload = json.load(f)

    job_id = payload["jobId"]
    files  = payload["files"]
    params = payload["params"]

    try:
        result = run(job_id, files, params)

        # The caller may hand us an explicit destination; otherwise fall back to
        # the portable temp directory (`/tmp` is not a valid path on Windows).
        result_path = args.output or os.path.join(
            tempfile.gettempdir(), f"result_{job_id}.json"
        )
        with open(result_path, "w", encoding="utf-8") as f:
            json.dump(result, f)

        log(f"Results written to {result_path}")
        sys.exit(0)
    except Exception as exc:  # noqa: BLE001 - top-level CLI guard: report and exit
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
