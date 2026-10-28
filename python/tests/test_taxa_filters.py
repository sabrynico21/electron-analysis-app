"""Real-data verification of the taxa filters (QMP fixture).

The app never requires a raw-count table: the correlation backends close the
data compositionally themselves and users may upload already-relative
abundances. The filters therefore have to evaluate *relative* abundances
regardless of the scale of the uploaded matrix. These tests pin that contract
against the committed QMP subset (20 samples x 30 taxa, SPRING::QMP):
prevalence, mean relative abundance and relative-abundance variance are compared
with an independent NumPy reference, and the removed raw-count filter is checked
to be inert even when a legacy config still carries it.

Run from the repository root::

    python -m pytest python/tests/test_taxa_filters.py -q
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from python import analysis

REPO_ROOT = Path(__file__).resolve().parents[2]
QMP_DIR = REPO_ROOT / 'spiec-easi-python' / 'tests' / 'golden' / 'qmp_n20_p30'


@pytest.fixture(scope='module')
def qmp():
    counts = np.loadtxt(QMP_DIR / 'counts.tsv', delimiter='\t')
    taxa = [line.strip() for line in (QMP_DIR / 'taxa.txt').read_text().splitlines() if line.strip()]
    samples = [line.strip() for line in (QMP_DIR / 'sample_ids.txt').read_text().splitlines() if line.strip()]

    assert counts.shape == (len(samples), len(taxa))
    assert np.all(counts >= 0)
    assert np.allclose(counts, np.round(counts)), 'QMP fixture is an integer-count matrix'

    reference = {
        # Presence is scale-free: a value above zero is the same after closure.
        'prevalence': (counts > 0).mean(axis=0),
        'relative': counts / counts.sum(axis=1, keepdims=True),
    }
    reference['mean_relative'] = reference['relative'].mean(axis=0)
    reference['variance'] = reference['relative'].var(axis=0)

    by_sample = {
        sample: {taxon: float(counts[row, column]) for column, taxon in enumerate(taxa)}
        for row, sample in enumerate(samples)
    }
    return {
        'counts': counts,
        'taxa': taxa,
        'samples': samples,
        'by_sample': by_sample,
        'reference': reference,
    }


def _config(**overrides):
    """Master toggle on, every individual filter off unless overridden."""
    filters = {
        'prevalence': {'enabled': False, 'percent': 20},
        'meanRelativeAbundance': {'enabled': False, 'percent': 0.01},
        'variability': {'enabled': False, 'minVariance': 0.01},
    }
    for name, values in overrides.items():
        filters[name].update(values)
    return {'taxaFiltering': {'enabled': True, 'preset': 'custom', **filters}}


def _apply(qmp, **overrides):
    restricted, kept, report = analysis._apply_taxa_filters(
        qmp['by_sample'],
        qmp['taxa'],
        'bacteria',
        _config(**overrides),
    )
    return restricted, kept, report


def _expected(taxa, mask):
    return [taxon for taxon, keep in zip(taxa, mask) if keep]


def test_qmp_prevalence_is_percentage_of_samples(qmp):
    _, kept, report = _apply(qmp, prevalence={'enabled': True, 'percent': 20})

    expected = _expected(qmp['taxa'], qmp['reference']['prevalence'] >= 0.20)
    assert kept == expected
    assert report['total_samples'] == len(qmp['samples']) == 20
    assert report['filters']['prevalence'] == {
        'enabled': True,
        'threshold_percent': 20.0,
        'removed': len(qmp['taxa']) - len(expected),
    }


def test_qmp_mean_relative_abundance_uses_percent_of_total_community(qmp):
    _, kept, report = _apply(qmp, meanRelativeAbundance={'enabled': True, 'percent': 0.01})

    mean_percent = qmp['reference']['mean_relative'] * 100.0
    expected = _expected(qmp['taxa'], mean_percent >= 0.01)
    assert kept == expected
    assert report['filters']['meanRelativeAbundance']['threshold_percent'] == pytest.approx(0.01)
    assert report['filters']['meanRelativeAbundance']['removed'] == len(qmp['taxa']) - len(expected)

    # The report must expose the same relative-abundance values the filter used.
    excluded = {row['taxon']: row for row in report['excluded']}
    for taxon, value in zip(qmp['taxa'], mean_percent):
        if taxon in excluded:
            assert excluded[taxon]['mean_relative_abundance_percent'] == pytest.approx(value, abs=5e-7)


def test_qmp_variability_is_variance_of_relative_abundance_not_raw_counts(qmp):
    threshold = 0.01
    _, kept, report = _apply(qmp, variability={'enabled': True, 'minVariance': threshold})

    reference_variance = qmp['reference']['variance']
    expected = _expected(qmp['taxa'], reference_variance >= threshold)
    assert kept == expected
    assert len(kept) == 3

    # If the filter had used the uploaded raw counts, every one of the 30 taxa
    # would pass this threshold; only 3 do on the relative scale.
    assert int((qmp['counts'].var(axis=0) >= threshold).sum()) == len(qmp['taxa'])
    assert report['filters']['variability']['threshold'] == pytest.approx(threshold)

    excluded = {row['taxon']: row for row in report['excluded']}
    for taxon, value in zip(qmp['taxa'], reference_variance):
        if taxon in excluded:
            assert excluded[taxon]['relative_abundance_variance'] == pytest.approx(value)


def test_qmp_combined_filters_are_and_ed(qmp):
    _, kept, report = _apply(
        qmp,
        prevalence={'enabled': True, 'percent': 20},
        meanRelativeAbundance={'enabled': True, 'percent': 0.01},
        variability={'enabled': True, 'minVariance': 0.0001},
    )

    reference = qmp['reference']
    expected = _expected(
        qmp['taxa'],
        (reference['prevalence'] >= 0.20)
        & (reference['mean_relative'] * 100.0 >= 0.01)
        & (reference['variance'] >= 0.0001),
    )
    assert kept == expected
    assert report['after'] == len(expected)
    # A taxon excluded by one filter is excluded even when it passes the others.
    for row in report['excluded']:
        assert set(row['reasons']) <= {'prevalence', 'meanRelativeAbundance', 'variability'}
        assert any(not keep for keep in (
            row['prevalence_percent'] >= 20.0,
            row['mean_relative_abundance_percent'] >= 0.01,
            row['relative_abundance_variance'] >= 0.0001,
        ))


def test_qmp_filters_are_scale_invariant_so_relative_input_is_accepted(qmp):
    """Raw counts, a per-sample rescale, and row-closed abundances agree."""
    combined = dict(
        prevalence={'enabled': True, 'percent': 20},
        meanRelativeAbundance={'enabled': True, 'percent': 0.01},
        variability={'enabled': True, 'minVariance': 0.0001},
    )
    _, kept_counts, report_counts = _apply(qmp, **combined)

    relative = qmp['reference']['relative']
    relative_by_sample = {
        sample: {taxon: float(relative[row, column]) for column, taxon in enumerate(qmp['taxa'])}
        for row, sample in enumerate(qmp['samples'])
    }
    _, kept_relative, report_relative = analysis._apply_taxa_filters(
        relative_by_sample, qmp['taxa'], 'bacteria', _config(**combined)
    )
    assert kept_relative == kept_counts

    # Rescaling each sample (different sequencing depth) must not change anything.
    depths = np.linspace(0.25, 40.0, len(qmp['samples']))[:, None]
    rescaled = qmp['counts'] * depths
    rescaled_by_sample = {
        sample: {taxon: float(rescaled[row, column]) for column, taxon in enumerate(qmp['taxa'])}
        for row, sample in enumerate(qmp['samples'])
    }
    _, kept_rescaled, _ = analysis._apply_taxa_filters(
        rescaled_by_sample, qmp['taxa'], 'bacteria', _config(**combined)
    )
    assert kept_rescaled == kept_counts

    # The reported relative metrics are identical across representations.
    for row_counts, row_relative in zip(report_counts['excluded'], report_relative['excluded']):
        assert row_counts['mean_relative_abundance_percent'] == pytest.approx(
            row_relative['mean_relative_abundance_percent'], abs=5e-7
        )
        assert row_counts['relative_abundance_variance'] == pytest.approx(
            row_relative['relative_abundance_variance'], rel=1e-9, abs=1e-12
        )


def test_removed_total_count_filter_is_inert(qmp):
    """A legacy `totalCount` entry must neither run nor appear in the report."""
    legacy = _config()
    legacy['taxaFiltering']['totalCount'] = {'enabled': True, 'minCount': 1e30}

    _, kept, report = analysis._apply_taxa_filters(
        qmp['by_sample'], qmp['taxa'], 'bacteria', legacy
    )

    assert kept == qmp['taxa']
    assert report['after'] == report['before'] == len(qmp['taxa'])
    assert 'totalCount' not in report['filters']
    assert all('total_count' not in row for row in report['excluded'])


def test_filtering_keeps_abundance_values_for_the_correlation_backends(qmp):
    """Filters select taxa; they do not replace the matrix with relative values."""
    restricted, kept, _ = _apply(qmp, prevalence={'enabled': True, 'percent': 20})

    for sample in qmp['samples']:
        assert set(restricted[sample]) == set(kept)
        for taxon in kept:
            assert restricted[sample][taxon] == qmp['by_sample'][sample][taxon]


# ── Report provenance: scope, denominator, guardrails ─────────────────────────


def test_report_declares_global_scope_domain_denominator_and_share(qmp):
    """The report states which samples and which denominator the filters used."""
    _, _, report = _apply(qmp, prevalence={'enabled': True, 'percent': 20})

    assert report['scope'] == 'all_samples_in_abundance_table'
    assert report['relative_abundance_denominator'] == 'domain'
    assert report['total_samples'] == len(qmp['samples'])
    expected_share = round(100.0 * (report['before'] - report['after']) / report['before'], 4)
    assert report['removed_share_percent'] == pytest.approx(expected_share, abs=1e-9)
    # The flag is derived from the unrounded share, so mirror that here.
    assert report['aggressive_removal'] is (
        (100.0 * (report['before'] - report['after']) / report['before']) > 60.0
    )


def test_disabled_filters_report_no_denominator_and_remove_nothing(qmp):
    _, kept, report = analysis._apply_taxa_filters(
        qmp['by_sample'], qmp['taxa'], 'bacteria', {'taxaFiltering': {'enabled': False}}
    )

    assert kept == qmp['taxa']
    assert report['enabled'] is False
    assert report['relative_abundance_denominator'] is None
    assert report['removed_share_percent'] == 0.0
    assert report['aggressive_removal'] is False
    assert report['scope'] == 'all_samples_in_abundance_table'


def test_aggressive_removal_is_flagged_for_the_user(qmp):
    """A threshold that wipes out the domain must be called out, not silent."""
    _, _, report = _apply(qmp, meanRelativeAbundance={'enabled': True, 'percent': 100})

    assert report['after'] == 0
    assert report['removed_share_percent'] == 100.0
    assert report['aggressive_removal'] is True


# ── Combined bacteria + fungi denominator (datasetType='both') ────────────────


def test_sample_totals_override_rescales_mean_and_variance_but_not_prevalence(qmp):
    """Doubling every per-sample total halves relative abundances.

    That is exactly the relationship between a per-domain denominator and the
    combined bacteria + fungi denominator the correlation model uses, so the
    thresholds move by 2x (mean) and 4x (variance) while prevalence — which is
    scale-free — must not move at all.
    """
    doubled_totals = {
        sample: 2.0 * float(counts_row.sum())
        for counts_row, sample in zip(qmp['counts'], qmp['samples'])
    }
    reference = qmp['reference']

    _, kept_mean, report_mean = analysis._apply_taxa_filters(
        qmp['by_sample'], qmp['taxa'], 'bacteria',
        _config(meanRelativeAbundance={'enabled': True, 'percent': 0.01}),
        sample_totals=doubled_totals,
    )
    # combined% = domain% / 2  ⇔  domain% ≥ 0.02
    assert kept_mean == _expected(qmp['taxa'], reference['mean_relative'] * 100.0 >= 0.02)
    assert report_mean['relative_abundance_denominator'] == 'combined_bacteria_fungi'
    excluded_mean = {row['taxon']: row for row in report_mean['excluded']}
    for taxon, value in zip(qmp['taxa'], reference['mean_relative'] * 50.0):
        if taxon in excluded_mean:
            assert excluded_mean[taxon]['mean_relative_abundance_percent'] == pytest.approx(
                value, abs=5e-7
            )

    _, kept_var, report_var = analysis._apply_taxa_filters(
        qmp['by_sample'], qmp['taxa'], 'bacteria',
        _config(variability={'enabled': True, 'minVariance': 0.04}),
        sample_totals=doubled_totals,
    )
    # variance scales with 1/4  ⇔  domain variance ≥ 0.16
    assert kept_var == _expected(qmp['taxa'], reference['variance'] >= 0.16)
    assert report_var['relative_abundance_denominator'] == 'combined_bacteria_fungi'

    _, kept_prevalence_default, _ = _apply(qmp, prevalence={'enabled': True, 'percent': 20})
    _, kept_prevalence_override, _ = analysis._apply_taxa_filters(
        qmp['by_sample'], qmp['taxa'], 'bacteria',
        _config(prevalence={'enabled': True, 'percent': 20}),
        sample_totals=doubled_totals,
    )
    assert kept_prevalence_override == kept_prevalence_default


def test_sample_totals_missing_entries_fall_back_to_domain_totals(qmp):
    """A sample the caller did not measure keeps its own domain denominator."""
    config = {'meanRelativeAbundance': {'enabled': True, 'percent': 0.01}}
    _, kept_default, report_default = _apply(qmp, **config)
    _, kept_fallback, report_fallback = analysis._apply_taxa_filters(
        qmp['by_sample'], qmp['taxa'], 'bacteria', _config(**config), sample_totals={}
    )

    assert kept_fallback == kept_default
    assert report_fallback['relative_abundance_denominator'] == 'combined_bacteria_fungi'
    for row_default, row_fallback in zip(report_default['excluded'], report_fallback['excluded']):
        assert row_fallback['mean_relative_abundance_percent'] == pytest.approx(
            row_default['mean_relative_abundance_percent'], abs=1e-9
        )


# ── Method-side removals (constant taxa / empty samples) ──────────────────────


def test_matrix_for_samples_records_method_side_drop_events():
    """Backend drops must reach the payload instead of living only in the log."""
    by_sample = {
        's1': {'varying': 1.0, 'flat': 2.0, 'zero': 0.0},
        's2': {'varying': 3.0, 'flat': 2.0, 'zero': 0.0},
        's3': {'varying': 5.0, 'flat': 2.0, 'zero': 0.0},
    }
    features = ['varying', 'flat', 'zero']
    samples = ['s1', 's2', 's3']

    events = []
    kept, matrix = analysis._matrix_for_samples(
        by_sample, features, samples, events, 'cclasso', 'bacteria · A'
    )

    assert kept == ['varying']
    assert matrix.shape == (3, 1)
    assert len(events) == 1
    event = events[0]
    assert event['method'] == 'cclasso'
    assert event['scope'] == 'bacteria · A'
    assert event['kind'] == 'constant_feature'
    assert event['item_kind'] == 'taxa'
    assert event['count'] == 2
    assert event['items'] == ['flat', 'zero']
    assert event['samples_in_scope'] == 3
    assert event['items_truncated'] is False

    # Without a collector the behaviour is unchanged (backwards compatible).
    kept_only, _ = analysis._matrix_for_samples(by_sample, features, samples)
    assert kept_only == ['varying']


def test_drop_event_payload_is_capped_and_optional_collector_is_a_no_op():
    events = []
    analysis._record_method_drop(
        events, 'sparcc', 'fungi · B', 'constant_feature',
        [f'taxon-{index}' for index in range(60)], 12,
    )

    [event] = events
    assert event['count'] == 60
    assert len(event['items']) == 50
    assert event['items_truncated'] is True
    assert event['item_kind'] == 'taxa'

    sample_events = []
    analysis._record_method_drop(
        sample_events, 'spring', 'bac+fun · A', 'empty_sample', ['s9'], 10
    )
    assert sample_events[0]['item_kind'] == 'samples'

    # No collector and no items ⇒ nothing is recorded.
    assert analysis._record_method_drop(None, 'sparcc', '', 'constant_feature', ['x'], 1) is None
    analysis._record_method_drop(events, 'sparcc', '', 'constant_feature', [], 1)
    assert len(events) == 1
