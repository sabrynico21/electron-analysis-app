import {
  DEFAULT_TAXA_FILTERING,
  TAXA_FILTER_PRESETS,
  TAXA_FILTER_PRESET_KEYS,
  applyTaxaFilteringPreset,
  describeTaxaFiltering,
  isTaxaFilteringValid,
  normalizeTaxaFiltering,
  prevalenceToMinSamples,
} from '../../utils/taxaFiltering'
import styles from './TaxaFilteringForm.module.css'

/** Editable numeric input with an inline unit and consistent typing semantics. */
function ThresholdInput({ id, value, onChange, onBlur, min, max, step, disabled, unit, 'aria-label': ariaLabel }) {
  return (
    <div className={styles.threshold}>
      <input
        id={id}
        type='number'
        inputMode='decimal'
        min={min}
        max={max}
        step={step}
        value={value}
        disabled={disabled}
        aria-label={ariaLabel}
        onChange={onChange}
        onBlur={onBlur}
      />
      {unit && <span className={styles.unit}>{unit}</span>}
    </div>
  )
}

/**
 * Taxa filtering configuration.
 *
 * Owns the `groupingMode` and `taxaFiltering` keys of the shared params object
 * (grouping was moved here from the correlation parameters because it belongs to
 * the taxa-preparation step). Every threshold is optional: a disabled master
 * toggle, or a disabled individual filter, means the corresponding threshold has
 * no effect on the run.
 */
export default function TaxaFilteringForm({
  values,
  onChange,
  disabled = false,
  totalSamples = 0,
}) {
  const config = normalizeTaxaFiltering(values.taxaFiltering)
  const clamp = (value, min, max) => Math.max(min, Math.min(max, value))

  const update = (patch) => {
    onChange((prev) => ({
      ...prev,
      taxaFiltering: { ...normalizeTaxaFiltering(prev.taxaFiltering), ...patch },
    }))
  }

  const updateFilter = (filterKey, patch) => {
    update({ [filterKey]: { ...config[filterKey], ...patch }, preset: 'custom' })
  }

  const handlePreset = (presetKey) => {
    update(applyTaxaFilteringPreset(config, presetKey))
  }

  const handleReset = () => {
    update(DEFAULT_TAXA_FILTERING)
  }

  const handleGroupingChange = (mode) => {
    onChange((prev) => ({ ...prev, groupingMode: mode }))
  }

  const masterDisabled = disabled || !config.enabled
  const filterDisabled = (filterKey) => masterDisabled || !config[filterKey].enabled
  const prevalenceMinSamples = prevalenceToMinSamples(config.prevalence.percent, totalSamples)
  const valid = isTaxaFilteringValid(values.taxaFiltering)
  // In combined mode the Python pipeline closes every sample on bacteria + fungi
  // together, so the mean/variability thresholds are measured on that same scale.
  const combinedDenominator = values.datasetType === 'both'
  const summary = describeTaxaFiltering(values.taxaFiltering, values.datasetType)

  /** Keep the input editable while typing; clamp/normalize only on blur. */
  const numberHandler = (filterKey, field, { min = 0, max = Infinity } = {}) => ({
    onChange: (e) => {
      const raw = e.target.value
      if (raw === '') {
        updateFilter(filterKey, { [field]: '' })
        return
      }
      const parsed = Number(raw)
      if (!Number.isFinite(parsed)) return
      const bounded = clamp(parsed, min, max)
      updateFilter(filterKey, { [field]: bounded })
    },
    onBlur: (e) => {
      const parsed = Number(e.target.value)
      const fallback = DEFAULT_TAXA_FILTERING[filterKey][field]
      const bounded = Number.isFinite(parsed) ? clamp(parsed, min, max) : fallback
      updateFilter(filterKey, { [field]: bounded })
    },
  })

  return (
    <div className={styles.container}>
      <div className={styles.masterRow}>
        <label className={styles.masterLabel} htmlFor='taxa-filtering-enabled'>
          <input
            id='taxa-filtering-enabled'
            type='checkbox'
            checked={config.enabled}
            disabled={disabled}
            onChange={(e) => update({ enabled: e.target.checked })}
          />
          <span>Enable taxa filtering</span>
        </label>
        <button
          type='button'
          className={styles.resetButton}
          onClick={handleReset}
          disabled={disabled}
          title='Restore the Balanced defaults'
        >
          Reset to defaults
        </button>
      </div>

      <div className={styles.field}>
        <label htmlFor='groupingMode'>Taxa grouping</label>
        <select
          id='groupingMode'
          value={values.groupingMode ?? 'all'}
          disabled={disabled}
          onChange={(e) => handleGroupingChange(e.target.value)}
        >
          <option value='all'>All taxa (no grouping)</option>
          <option value='phylum'>Group by phylum (p_)</option>
          <option value='genus'>Group by genus (g_)</option>
        </select>
        <small className={styles.help}>Grouping is applied after filtering, on the surviving taxa.</small>
      </div>

      <div className={styles.presetsBlock}>
        <span className={styles.presetsLabel}>Presets</span>
        <div className={styles.presets} role='group' aria-label='Taxa filtering presets'>
          {TAXA_FILTER_PRESET_KEYS.map((key) => {
            const preset = TAXA_FILTER_PRESETS[key]
            const active = config.preset === key
            return (
              <button
                key={key}
                type='button'
                className={`${styles.presetButton} ${active ? styles.presetActive : ''}`}
                onClick={() => handlePreset(key)}
                disabled={disabled}
                aria-pressed={active}
                title={preset.description}
              >
                {preset.label}
              </button>
            )
          })}
        </div>
        <small className={styles.help}>
          Presets are starting points, not universally optimal values.
        </small>
      </div>

      <div className={styles.filters}>
        {/* A. Minimum prevalence */}
        <div className={`${styles.filterCard} ${!config.enabled ? styles.filterInactive : ''}`}>
          <label className={styles.filterHeader} htmlFor='filter-prevalence-enabled'>
            <input
              id='filter-prevalence-enabled'
              type='checkbox'
              checked={config.prevalence.enabled}
              disabled={disabled || !config.enabled}
              onChange={(e) => updateFilter('prevalence', { enabled: e.target.checked })}
            />
            <span className={styles.filterTitle}>Minimum prevalence</span>
          </label>
          <p className={styles.filterHelp}>
            Keeps taxa present (abundance &gt; 0) in at least this percentage of samples.
            Evaluated on every sample of the abundance table.
          </p>
          <ThresholdInput
            id='filter-prevalence-percent'
            aria-label='Minimum prevalence percentage'
            min={0}
            max={100}
            step={1}
            unit='%'
            value={config.prevalence.percent}
            disabled={filterDisabled('prevalence')}
            {...numberHandler('prevalence', 'percent', { min: 0, max: 100 })}
          />
          {totalSamples > 0 && (
            <small className={styles.derived}>
              ≈ {prevalenceMinSamples} of {totalSamples} samples (metadata count)
            </small>
          )}
        </div>

        {/* B. Minimum mean relative abundance */}
        <div className={`${styles.filterCard} ${!config.enabled ? styles.filterInactive : ''}`}>
          <label className={styles.filterHeader} htmlFor='filter-meanra-enabled'>
            <input
              id='filter-meanra-enabled'
              type='checkbox'
              checked={config.meanRelativeAbundance.enabled}
              disabled={disabled || !config.enabled}
              onChange={(e) => updateFilter('meanRelativeAbundance', { enabled: e.target.checked })}
            />
            <span className={styles.filterTitle}>Minimum mean relative abundance</span>
          </label>
          <p className={styles.filterHelp}>
            Keeps taxa whose mean relative abundance across all samples is at least this value.
            Absent samples count as 0.
            {combinedDenominator
              ? ' Relative to the combined bacteria + fungi total — the same scale the correlation model uses.'
              : ' Relative to the total of this domain.'}
          </p>
          <ThresholdInput
            id='filter-meanra-percent'
            aria-label='Minimum mean relative abundance percentage'
            min={0}
            step={0.001}
            unit='%'
            value={config.meanRelativeAbundance.percent}
            disabled={filterDisabled('meanRelativeAbundance')}
            {...numberHandler('meanRelativeAbundance', 'percent', { min: 0 })}
          />
        </div>

        {/* C. Minimum variability */}
        <div className={`${styles.filterCard} ${!config.enabled ? styles.filterInactive : ''}`}>
          <label className={styles.filterHeader} htmlFor='filter-variability-enabled'>
            <input
              id='filter-variability-enabled'
              type='checkbox'
              checked={config.variability.enabled}
              disabled={disabled || !config.enabled}
              onChange={(e) => updateFilter('variability', { enabled: e.target.checked })}
            />
            <span className={styles.filterTitle}>Minimum variability</span>
          </label>
          <p className={styles.filterHelp}>
            Keeps taxa whose relative abundance varies across samples (variance ≥ threshold). Removes near-constant taxa.
            {combinedDenominator
              ? ' Variance is measured on the combined bacteria + fungi scale.'
              : ''}
          </p>
          <ThresholdInput
            id='filter-variability-min'
            aria-label='Minimum relative-abundance variance'
            min={0}
            step={0.0001}
            unit='var'
            value={config.variability.minVariance}
            disabled={filterDisabled('variability')}
            {...numberHandler('variability', 'minVariance', { min: 0 })}
          />
        </div>
      </div>

      <div className={styles.summary} aria-live='polite'>
        <span className={styles.summaryTitle}>Active filters</span>
        {!config.enabled ? (
          <p className={styles.summaryOff}>
            Taxa filtering is off — no user-defined threshold is applied. Constant and all-zero
            taxa are still discarded by the correlation method.
          </p>
        ) : (
          <ul className={styles.summaryList}>
            {summary.lines.map((line) => (
              <li key={line}>{line}</li>
            ))}
          </ul>
        )}
      </div>

      <p className={styles.help}>
        The correlation method always discards constant and
        all-zero taxa (and SPRING also discards empty samples).
      </p>

      {!valid && (
        <p className={styles.validation} role='alert'>
          Some thresholds are invalid. Use 0–100% for prevalence, non-negative values elsewhere.
        </p>
      )}

    </div>
  )
}
