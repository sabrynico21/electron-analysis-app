/**
 * Deterministic fixtures for the end-to-end suites.
 *
 * Every suite runs against an isolated Electron user-data directory generated
 * here, so the tests are hermetic: they never read or write the real archive and
 * they do not need any prior analysis to have been run.
 *
 * Layout produced (mirrors what the application itself writes):
 *   <userData>/analysis-history.json
 *   <userData>/analysis-results/<jobId>.json
 *   <userData>/analysis-cache/index.json
 *   <userData>/analysis-cache/graphs/<cacheKey>.json
 *   <inputDir>/{abundance.csv,metadata.tsv,metadata_no_conditions.tsv,metadata_no_column.tsv}
 */
import fs from 'node:fs'
import path from 'node:path'

const CACHE_VERSION = 1

/** Fixed timestamps keep the history order and the rendered dates deterministic. */
export const TIMESTAMPS = {
  richAnalysis: Date.parse('2026-01-05T10:00:00Z'),
  plainAnalysis: Date.parse('2026-01-04T09:00:00Z'),
  missingFilesAnalysis: Date.parse('2026-01-03T08:00:00Z'),
}

export const CONDITION_A = 'Control'
export const CONDITION_B = 'Treated'

export const IDS = {
  rich: 'aaaaaaaa-0000-4000-8000-000000000001',
  plain: 'bbbbbbbb-0000-4000-8000-000000000002',
  missingFiles: 'cccccccc-0000-4000-8000-000000000003',
}

const FEATURES_CONDITION_A = ['otu_alpha', 'otu_beta', 'otu_gamma']
const FEATURES_CONDITION_B = ['otu_beta', 'otu_delta']

/** Input tables written to disk so the upload pipeline has real files to read. */
function writeInputFiles(inputDir) {
  fs.mkdirSync(inputDir, { recursive: true })

  fs.writeFileSync(
    path.join(inputDir, 'abundance.csv'),
    [
      'sample-id,otu_alpha,otu_beta,otu_gamma,otu_delta',
      'S1,12,4,31,7',
      'S2,15,6,28,9',
      'S3,9,3,35,6',
      'S4,18,8,24,11',
      'S5,11,5,33,8',
      'S6,14,7,26,10',
      'S7,21,2,38,4',
      'S8,24,9,19,13',
      '',
    ].join('\n'),
    'utf-8'
  )

  const metadataRows = (conditions) => [
    ['sample-id', 'Combination_treat2', 'Group'],
    ...conditions.map((condition, index) => [`S${index + 1}`, condition, index < 4 ? '1' : '2']),
  ]
  const toTsv = (rows) => `${rows.map((row) => row.join('\t')).join('\n')}\n`

  fs.writeFileSync(
    path.join(inputDir, 'metadata.tsv'),
    toTsv(metadataRows([
      CONDITION_A, CONDITION_A, CONDITION_A, CONDITION_A,
      CONDITION_B, CONDITION_B, CONDITION_B, CONDITION_B,
    ])),
    'utf-8'
  )

  fs.writeFileSync(
    path.join(inputDir, 'metadata_semicolon.csv'),
    ['sample-id;Combination_treat2', 'S1;A', 'S2;B'].join('\n') + '\n',
    'utf-8'
  )

  // Only one distinct value in Combination_treat2 -> must be rejected.
  fs.writeFileSync(
    path.join(inputDir, 'metadata_no_conditions.tsv'),
    toTsv([
      ['sample-id', 'Combination_treat2'],
      ['S1', 'OnlyOne'],
      ['S2', 'OnlyOne'],
    ]),
    'utf-8'
  )

  // Missing the required column -> must be rejected with an explicit message.
  fs.writeFileSync(
    path.join(inputDir, 'metadata_no_column.tsv'),
    toTsv([
      ['sample-id', 'SomeOtherColumn'],
      ['S1', 'A'],
      ['S2', 'B'],
    ]),
    'utf-8'
  )

  // Same content as abundance.csv but with a different name, used to check that
  // re-adding the same path is de-duplicated.
  fs.writeFileSync(
    path.join(inputDir, 'abundance_copy.csv'),
    fs.readFileSync(path.join(inputDir, 'abundance.csv'), 'utf-8'),
    'utf-8'
  )
}

/** A condition graph with a handful of nodes and edges. */
function conditionGraph(features) {
  const nodes = features.map((id) => ({ id, label: id, kind: 'bacteria' }))
  const edges = features.slice(1).map((target, index) => ({
    source: features[index],
    target,
    weight: index % 2 === 0 ? 0.6 : -0.4,
    p_value: 0.01,
    method: 'spearman',
  }))
  return {
    graph_id: 'graph-fixture',
    relation_type: 'condition',
    nodes,
    edges,
    stats: { node_count: nodes.length, edge_count: edges.length },
  }
}

/** A saved clustering result, shaped like the response of run_clustering.py. */
function clustering(seedNode, features) {
  const nodes = features.map((id) => ({ id, label: id, kind: 'bacteria' }))
  const edges = features.slice(1).map((target, index) => ({
    source: features[index],
    target,
    weight: 0.5,
    p_value: 0.02,
    method: 'spearman',
  }))
  return {
    success: true,
    cacheKey: 'placeholder',
    seedNode,
    condition: 'placeholder',
    algorithm: 'pagerank_nibble_adaptive',
    clusters: [features],
    clusterNodes: features,
    subgraph: { nodes, edges },
    metadata: {
      b: 3,
      phi_target: 0.2,
      p_at_vol: 0.42,
      weighted_mode: true,
      ranked_nodes: features.map((id, index) => ({ id, rank: index })),
    },
    message: 'fixture',
  }
}

/**
 * Build one analysis: a result file, its cached condition graphs and the cache
 * index entries. Returns the history payload snapshot and the graph refs.
 */
function writeAnalysis({ userDataDir, jobId, name, createdAt, files, withClusters, cacheIndex }) {
  const resultsDir = path.join(userDataDir, 'analysis-results')
  const graphsDir = path.join(userDataDir, 'analysis-cache', 'graphs')
  fs.mkdirSync(resultsDir, { recursive: true })
  fs.mkdirSync(graphsDir, { recursive: true })

  const conditions = [
    { key: 'cond-a', condition: CONDITION_A, features: FEATURES_CONDITION_A, seed: 'otu_beta' },
    { key: 'cond-b', condition: CONDITION_B, features: FEATURES_CONDITION_B, seed: 'otu_delta' },
  ]

  const graphRefs = []
  const lastClustering = {}

  conditions.forEach(({ key, condition, features, seed }) => {
    const cacheKey = `${jobId.slice(0, 8)}-${key}`
    const graph = conditionGraph(features)
    fs.writeFileSync(
      path.join(graphsDir, `${cacheKey}.json`),
      JSON.stringify({
        version: CACHE_VERSION,
        cacheKey,
        jobId,
        condition,
        correlationMethod: 'spearman',
        datasetType: 'bacteria',
        pValueThreshold: 0.05,
        createdAt: new Date(createdAt).toISOString(),
        graph,
        pairFile: null,
      }),
      'utf-8'
    )
    const ref = { cacheKey, path: path.join(graphsDir, `${cacheKey}.json`), version: CACHE_VERSION, jobId, condition }
    graphRefs.push(ref)
    cacheIndex[cacheKey] = ref

    if (withClusters) {
      lastClustering[cacheKey] = { ...clustering(seed, features), cacheKey, condition }
    }
  })

  const summary = {
    n_bacteria_files: (files.bacteria || []).length,
    n_bacteria_metadata_files: (files.bacteria_metadata || []).length,
    dataset_type: 'bacteria',
    correlation_method: 'spearman',
    condition_a: CONDITION_A,
    condition_b: CONDITION_B,
    grouping_mode: 'all',
    p_value_threshold: 0.05,
    p_value: 0.0312,
    statistic: withClusters ? 4.721 : 2.5,
    status: 'significant',
  }

  const result = {
    summary,
    charts: [],
    table: [
      { row: 'otu_alpha', column: 'otu_beta', cor: 0.62, p: 0.01, method: 'spearman', condition: CONDITION_A },
      { row: 'otu_beta', column: 'otu_gamma', cor: -0.41, p: 0.02, method: 'spearman', condition: CONDITION_A },
      { row: 'otu_beta', column: 'otu_delta', cor: 0.33, p: 0.03, method: 'spearman', condition: CONDITION_B },
    ],
    graph_pipeline: {
      version: 1,
      status: 'ok',
      condition_graphs: conditions.map(({ condition, features }) => ({
        graph_id: 'graph-fixture',
        condition,
        correlation_method: 'spearman',
        pair_file: { format: 'tsv', path: null, row_count: features.length },
        graph: conditionGraph(features),
        clustering: { status: 'not_run', algorithm: null, seed_node: null, clusters: [], message: '' },
      })),
      cache_version: CACHE_VERSION,
      graph_refs: graphRefs,
      last_clustering: lastClustering,
    },
    input: {
      jobId,
      scriptName: 'analysis.py',
      analysisName: name,
      files,
      params: {
        datasetType: 'bacteria',
        correlationMethod: 'spearman',
        groupingMode: 'all',
        pValueThreshold: 0.05,
        conditionA: CONDITION_A,
        conditionB: CONDITION_B,
      },
      completedAt: createdAt,
    },
  }

  fs.writeFileSync(path.join(resultsDir, `${jobId}.json`), JSON.stringify(result), 'utf-8')

  return {
    entry: {
      id: jobId,
      analysisName: name,
      scriptName: 'analysis.py',
      status: 'completed',
      createdAt,
      updatedAt: createdAt,
      runCount: 1,
      payloadSnapshot: { files, params: result.input.params },
      summary,
      hasResults: true,
    },
  }
}

/**
 * Create an isolated user-data directory plus the input files the tests upload.
 * @param {string} baseDir a temporary directory owned by the caller
 */
export function createFixtures(baseDir) {
  const userDataDir = path.join(baseDir, 'user-data')
  const inputDir = path.join(baseDir, 'inputs')
  fs.mkdirSync(userDataDir, { recursive: true })
  writeInputFiles(inputDir)

  const files = {
    bacteria: [path.join(inputDir, 'abundance.csv')],
    bacteria_metadata: [path.join(inputDir, 'metadata.tsv')],
    fungi: [],
    fungi_metadata: [],
  }

  const cacheIndex = {}
  const rich = writeAnalysis({
    userDataDir, jobId: IDS.rich, name: 'Rich analysis', createdAt: TIMESTAMPS.richAnalysis, files, withClusters: true, cacheIndex,
  })
  const plain = writeAnalysis({
    userDataDir, jobId: IDS.plain, name: 'Plain analysis', createdAt: TIMESTAMPS.plainAnalysis, files, withClusters: false, cacheIndex,
  })
  const missingFiles = writeAnalysis({
    userDataDir, jobId: IDS.missingFiles, name: 'Moved files analysis', createdAt: TIMESTAMPS.missingFilesAnalysis, files: { bacteria: [path.join(inputDir, 'gone.csv')] }, withClusters: false, cacheIndex,
  })

  fs.mkdirSync(path.join(userDataDir, 'analysis-cache'), { recursive: true })
  fs.writeFileSync(
    path.join(userDataDir, 'analysis-cache', 'index.json'),
    JSON.stringify(cacheIndex, null, 2),
    'utf-8'
  )

  fs.writeFileSync(
    path.join(userDataDir, 'analysis-history.json'),
    JSON.stringify({
      version: 2,
      updatedAt: new Date().toISOString(),
      entries: [rich.entry, plain.entry, missingFiles.entry],
    }, null, 2),
    'utf-8'
  )

  return { userDataDir, inputDir, files, entries: [rich.entry, plain.entry, missingFiles.entry] }
}
