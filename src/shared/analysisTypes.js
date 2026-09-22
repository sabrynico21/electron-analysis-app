/**
 * Shared type definitions (JSDoc) used by both main and renderer.
 *
 * @typedef {'idle'|'running'|'completed'|'failed'|'cancelled'} JobStatus
 *
 * @typedef {Object} AnalysisPayload
 * @property {string}   scriptName
 * @property {string[]} files
 * @property {Object}   params
 *
 * @typedef {Object} ProgressEvent
 * @property {string}    jobId
 * @property {JobStatus} status
 * @property {number}    percent
 *
 * @typedef {Object} AnalysisResult
 * @property {Object}   summary   - Key-value stats
 * @property {Object[]} charts    - Chart.js datasets
 * @property {Object[]} table     - Tabular rows
 * @property {GraphPipeline=} graph_pipeline - Condition graph build artifacts
 *
 * @typedef {Object} GraphNode
 * @property {string} id
 * @property {string} label
 * @property {'bacteria'|'fungi'} kind
 *
 * @typedef {Object} GraphEdge
 * @property {string} source
 * @property {string} target
 * @property {number} weight
 * @property {number} p_value
 * @property {string} method
 *
 * @typedef {Object} ConditionGraph
 * @property {string} graph_id
 * @property {string} condition
 * @property {string} correlation_method
 * @property {{ format: string, path: string, row_count: number }} pair_file
 * @property {{ nodes: GraphNode[], edges: GraphEdge[], stats: Object }} graph
 * @property {{ status: string, algorithm: (string|null), seed_node: (string|null), clusters: Object[], message: string }} clustering
 *
 * @typedef {Object} GraphPipeline
 * @property {number} version
 * @property {string} status
 * @property {ConditionGraph[]} condition_graphs
 */
