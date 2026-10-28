/**
 * Taxa filtering configuration: defaults, presets, validation and helpers.
 *
 * The whole configuration lives under `params.taxaFiltering` so it travels with
 * the analysis payload (and the history snapshot) untouched, exactly like the
 * other analysis parameters. Keeping the shape in one place means the renderer,
 * the history restore path and the Python pipeline all agree on the contract.
 *
 * @typedef {Object} TaxaFilter
 * @property {boolean} enabled
 * @property {number|string} value  Threshold; a string while the user is typing.
 *
 * @typedef {Object} TaxaFilteringConfig
 * @property {boolean} enabled
 * @property {string}  preset
 * @property {{enabled: boolean, percent: number|string}} prevalence
 * @property {{enabled: boolean, percent: number|string}} meanRelativeAbundance
 * @property {{enabled: boolean, minVariance: number|string}} variability
 */

/** Filtering off altogether — no taxa are removed. */
export const DEFAULT_TAXA_FILTERING = {
  enabled: true,
  preset: 'balanced',
  // A. Minimum prevalence, as a percentage of samples with abundance > 0.
  prevalence: { enabled: true, percent: 20 },
  // B. Minimum mean relative abundance across all samples, as a percentage.
  meanRelativeAbundance: { enabled: true, percent: 0.01 },
  // C. Minimum variability: variance of relative abundances across samples.
  variability: { enabled: false, minVariance: 0.01 },
}

/**
 * Editable presets. They are starting points, not universally optimal values:
 * every threshold stays editable after a preset is applied.
 */
export const TAXA_FILTER_PRESETS = {
  exploratory: {
    label: 'Exploratory',
    description: 'Keeps more taxa — good for an initial look at the data.',
    config: {
      prevalence: { enabled: true, percent: 10 },
      meanRelativeAbundance: { enabled: false, percent: 0.01 },
      variability: { enabled: false },
    },
  },
  balanced: {
    label: 'Balanced',
    description: 'Recommended default — removes rare taxa without being aggressive.',
    config: {
      prevalence: { enabled: true, percent: 20 },
      meanRelativeAbundance: { enabled: true, percent: 0.01 },
      variability: { enabled: false },
    },
  },
  stringent: {
    label: 'Stringent',
    description: 'Keeps only well-represented taxa — fewer, more robust features.',
    config: {
      prevalence: { enabled: true, percent: 50 },
      meanRelativeAbundance: { enabled: true, percent: 0.05 },
      variability: { enabled: false },
    },
  },
}

const isPlainObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value)

/**
 * Build a complete configuration from a (possibly partial/legacy) one.
 *
 * Presets intentionally keep the numeric value of the filters they disable so a
 * user who flips a toggle back on finds their previous threshold intact.
 */
export function normalizeTaxaFiltering(raw) {
  const source = isPlainObject(raw) ? raw : {}
  const prevalence = isPlainObject(source.prevalence) ? source.prevalence : {}
  const meanRelativeAbundance = isPlainObject(source.meanRelativeAbundance) ? source.meanRelativeAbundance : {}
  const variability = isPlainObject(source.variability) ? source.variability : {}

  const num = (value, fallback) => (value === '' || value === null || value === undefined ? fallback : Number(value))

  return {
    enabled: source.enabled ?? DEFAULT_TAXA_FILTERING.enabled,
    preset: typeof source.preset === 'string' ? source.preset : DEFAULT_TAXA_FILTERING.preset,
    prevalence: {
      enabled: prevalence.enabled ?? DEFAULT_TAXA_FILTERING.prevalence.enabled,
      percent: num(prevalence.percent, DEFAULT_TAXA_FILTERING.prevalence.percent),
    },
    meanRelativeAbundance: {
      enabled: meanRelativeAbundance.enabled ?? DEFAULT_TAXA_FILTERING.meanRelativeAbundance.enabled,
      percent: num(meanRelativeAbundance.percent, DEFAULT_TAXA_FILTERING.meanRelativeAbundance.percent),
    },
    variability: {
      enabled: variability.enabled ?? DEFAULT_TAXA_FILTERING.variability.enabled,
      minVariance: num(variability.minVariance, DEFAULT_TAXA_FILTERING.variability.minVariance),
    },
  }
}

/**
 * Apply a preset on top of the current configuration.
 *
 * Only the filters the preset mentions are updated; the rest (and every numeric
 * value the preset disables) is preserved so switching preset never destroys a
 * custom value the user may want to re-enable.
 */
export function applyTaxaFilteringPreset(current, presetKey) {
  const preset = TAXA_FILTER_PRESETS[presetKey]
  const base = normalizeTaxaFiltering(current)
  if (!preset) return base

  const next = { ...base, enabled: true, preset: presetKey }
  for (const [filterKey, override] of Object.entries(preset.config)) {
    next[filterKey] = { ...base[filterKey], ...override }
  }
  return next
}

/** The presets as an array, for rendering. */
export const TAXA_FILTER_PRESET_KEYS = ['exploratory', 'balanced', 'stringent']

/**
 * Whether the configuration is safe to run.
 *
 * A disabled master toggle (or a disabled individual filter) makes its threshold
 * irrelevant, so it is never validated in that state.
 */
export function isTaxaFilteringValid(raw) {
  const config = normalizeTaxaFiltering(raw)
  if (!config.enabled) return true

  if (config.prevalence.enabled) {
    const p = Number(config.prevalence.percent)
    if (!Number.isFinite(p) || p < 0 || p > 100) return false
  }
  if (config.meanRelativeAbundance.enabled) {
    const p = Number(config.meanRelativeAbundance.percent)
    if (!Number.isFinite(p) || p < 0) return false
  }
  if (config.variability.enabled) {
    const v = Number(config.variability.minVariance)
    if (!Number.isFinite(v) || v < 0) return false
  }
  return true
}

/**
 * Human-readable summary of the active settings, used for the live readout.
 *
 * The wording is deliberately precise about two things the numbers alone hide:
 *
 * - "off" never means "every taxon is used": the correlation backends always
 *   discard constant / all-zero taxa (zero variance cannot produce a
 *   correlation), so we say so instead of promising the full taxon set.
 * - In combined mode (`datasetType === 'both'`) the mean/variability thresholds
 *   are evaluated on the combined bacteria + fungi per-sample total, which is
 *   the scale the correlation model actually uses (see `_apply_taxa_filters`).
 *
 * @param {Object}  raw         Current (possibly partial) configuration.
 * @param {string}  [datasetType]  'bacteria' | 'fungi' | 'both'.
 * @returns {{active: boolean, lines: string[]}}
 */
export function describeTaxaFiltering(raw, datasetType) {
  const config = normalizeTaxaFiltering(raw)
  const combinedDenominator = datasetType === 'both'
  const constantTaxaNote =
    'Constant and all-zero taxa are still discarded by the correlation method — they cannot produce a correlation.'

  if (!config.enabled) {
    return {
      active: false,
      lines: ['Taxa filtering is off — no user-defined threshold is applied.', constantTaxaNote],
    }
  }

  const lines = []
  if (config.prevalence.enabled) {
    lines.push(`Prevalence ≥ ${formatNumber(config.prevalence.percent)}% of samples`)
  }
  if (config.meanRelativeAbundance.enabled) {
    lines.push(
      combinedDenominator
        ? `Mean relative abundance ≥ ${formatNumber(config.meanRelativeAbundance.percent)}% of the combined bacteria + fungi total`
        : `Mean relative abundance ≥ ${formatNumber(config.meanRelativeAbundance.percent)}%`
    )
  }
  if (config.variability.enabled) {
    lines.push(
      combinedDenominator
        ? `Relative-abundance variance ≥ ${formatNumber(config.variability.minVariance)} (combined bacteria + fungi scale)`
        : `Relative-abundance variance ≥ ${formatNumber(config.variability.minVariance)}`
    )
  }

  const noIndividualFilter = lines.length === 0
  if (noIndividualFilter) {
    return {
      active: false,
      lines: [
        'No individual filter is active — no user-defined threshold is applied.',
        constantTaxaNote,
      ],
    }
  }

  // Scope is global on purpose: it keeps the candidate taxon set identical for
  // both selected conditions (and therefore the same graph node pool).
  lines.push('Scope: every sample of the abundance table (all conditions).')
  return { active: true, lines }
}

/** Compact number formatting that never shows trailing noise for small values. */
export function formatNumber(value) {
  const n = Number(value)
  if (!Number.isFinite(n)) return String(value ?? '')
  if (n === 0) return '0'
  if (Math.abs(n) < 0.01) return n.toPrecision(1)
  return String(Number(n.toFixed(4)))
}

/** Equivalent minimum sample count for a prevalence percentage. */
export function prevalenceToMinSamples(percent, totalSamples) {
  const p = Number(percent)
  const total = Number(totalSamples)
  if (!Number.isFinite(p) || !Number.isFinite(total) || total <= 0) return null
  return Math.ceil((p / 100) * total)
}
