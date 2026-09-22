export const GROUPING_MODE = {
  ALL: 'all',
  PHYLUM: 'phylum',
  GENUS: 'genus',
}

function tokenizeTaxonomy(value) {
  return String(value || '')
    .split(/[|;,\s]+/)
    .map((part) => part.trim())
    .filter(Boolean)
}

function findTokenByPrefix(value, prefix) {
  const tokens = tokenizeTaxonomy(value)
  return tokens.find((token) => token.toLowerCase().startsWith(prefix)) || ''
}

function normalizeToken(token, fallback) {
  if (!token) return fallback
  const cleaned = token.replace(/^[a-z]__?/i, '').trim()
  return cleaned || fallback
}

export function getNodeDisplayName(node, groupingMode = GROUPING_MODE.ALL) {
  const raw = node?.label || node?.id || ''

  if (groupingMode === GROUPING_MODE.PHYLUM) {
    const token = findTokenByPrefix(raw, 'p_')
    return normalizeToken(token, 'phylum_unknown')
  }

  if (groupingMode === GROUPING_MODE.GENUS) {
    const token = findTokenByPrefix(raw, 'g_')
    return normalizeToken(token, 'genus_unknown')
  }

  return String(raw || node?.id || 'otu_unknown')
}

export function getGroupingNodeKey(node, groupingMode = GROUPING_MODE.ALL) {
  if (groupingMode === GROUPING_MODE.PHYLUM) {
    return `phylum:${getNodeDisplayName(node, groupingMode)}`
  }
  if (groupingMode === GROUPING_MODE.GENUS) {
    return `genus:${getNodeDisplayName(node, groupingMode)}`
  }
  return String(node?.id || '')
}
