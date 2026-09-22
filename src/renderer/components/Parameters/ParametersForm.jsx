import styles from './ParametersForm.module.css'

export default function ParametersForm({ values, onChange, availableDatasets, conditionOptions, disabled = false }) {
  const handleChange = (key, value) => onChange((prev) => ({ ...prev, [key]: value }))
  const clamp = (value, min, max) => Math.max(min, Math.min(max, value))

  const datasetOptions = [
    { value: 'bacteria', label: 'Only bacteria', enabled: availableDatasets.bacteria },
    { value: 'fungi', label: 'Only fungi', enabled: availableDatasets.fungi },
    { value: 'both', label: 'Bacteria + fungi', enabled: availableDatasets.bacteria && availableDatasets.fungi },
  ]

  return (
    <form className={styles.form} onSubmit={(e) => e.preventDefault()}>
      <div className={styles.field}>
        <label htmlFor="datasetType">Data type to analyze</label>
        <select
          id="datasetType"
          value={values.datasetType ?? ''}
          disabled={disabled}
          onChange={(e) => handleChange('datasetType', e.target.value)}
        >
          <option value="" disabled>Select available data type</option>
          {datasetOptions.map((opt) => (
            <option key={opt.value} value={opt.value} disabled={!opt.enabled}>{opt.label}</option>
          ))}
        </select>
      </div>

      <div className={styles.field}>
        <label htmlFor="correlationMethod">Correlation method</label>
        <select
          id="correlationMethod"
          value={values.correlationMethod ?? 'spearman'}
          disabled={disabled}
          onChange={(e) => handleChange('correlationMethod', e.target.value)}
        >
          <option value="spearman">Spearman</option>
          <option value="sparcc">SparCC</option>
        </select>
      </div>

      <div className={styles.field}>
        <label htmlFor="groupingMode">Taxa grouping</label>
        <select
          id="groupingMode"
          value={values.groupingMode ?? 'all'}
          disabled={disabled}
          onChange={(e) => handleChange('groupingMode', e.target.value)}
        >
          <option value="all">All taxa (no grouping)</option>
          <option value="phylum">Group by phylum (p_)</option>
          <option value="genus">Group by genus (g_)</option>
        </select>
      </div>

      <div className={styles.field}>
        <label htmlFor="pValueThreshold">
          P-value threshold <span className={styles.range}>(0.0 - 1.0)</span>
        </label>
        <input
          id="pValueThreshold"
          type="number"
          min="0"
          max="1"
          step="0.001"
          disabled={disabled}
          value={values.pValueThreshold ?? 0.05}
          onChange={(e) => {
            // Mantieni input editabile ma non permettere di superare i limiti.
            const raw = e.target.value
            if (raw === '') {
              handleChange('pValueThreshold', '')
              return
            }
            const parsed = Number(raw)
            if (Number.isFinite(parsed)) {
              handleChange('pValueThreshold', clamp(parsed, 0, 1))
            }
          }}
          onBlur={(e) => {
            // Valida e applica i limiti solo quando l'utente smette di digitare
            const parsed = Number(e.target.value)
            const normalized = Number.isFinite(parsed) ? clamp(parsed, 0, 1) : 0.05
            handleChange('pValueThreshold', normalized)
          }}
        />
      </div>

      <div className={styles.field}>
        <label htmlFor="conditionA">Condition A</label>
        <select
          id="conditionA"
          value={values.conditionA ?? ''}
          onChange={(e) => handleChange('conditionA', e.target.value)}
          disabled={disabled || conditionOptions.length === 0}
        >
          <option value="" disabled>Select first condition</option>
          {conditionOptions.map((opt) => (
            <option key={opt} value={opt}>{opt}</option>
          ))}
        </select>
      </div>

      <div className={styles.field}>
        <label htmlFor="conditionB">Condition B</label>
        <select
          id="conditionB"
          value={values.conditionB ?? ''}
          onChange={(e) => handleChange('conditionB', e.target.value)}
          disabled={disabled || conditionOptions.length === 0}
        >
          <option value="" disabled>Select second condition</option>
          {conditionOptions.map((opt) => (
            <option key={opt} value={opt} disabled={opt === values.conditionA}>{opt}</option>
          ))}
        </select>
      </div>
    </form>
  )
}
