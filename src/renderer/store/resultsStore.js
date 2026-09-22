import { create } from 'zustand'

export const useResultsStore = create((set) => ({
  results: null,
  conditionGraphs: {},
  reclusterByGraph: {},
  isLoading: false,
  error: null,

  fetchResults: async (jobId) => {
    set({ isLoading: true, error: null })
    try {
      const data = await window.electronAPI.getResults(jobId)
      if (data) {
        set({ results: data, isLoading: false, conditionGraphs: {}, reclusterByGraph: {} })
        return data
      }
      set({ results: null, error: 'Results not found for this analysis', isLoading: false })
      return null
    } catch (error) {
      set({ results: null, error: error?.message || 'Unable to load results', isLoading: false })
      return null
    }
  },

  fetchConditionGraph: async (cacheKey) => {
    if (!cacheKey) return
    const response = await window.electronAPI.getConditionGraph(cacheKey)
    if (!response?.success) {
      set({ error: response?.error || 'Unable to load condition graph' })
      return
    }

    set((state) => ({
      conditionGraphs: {
        ...state.conditionGraphs,
        [cacheKey]: response.graph,
      },
    }))
  },

  reclusterConditionGraph: async ({ cacheKey, seedNode, options = {} }) => {
    if (!cacheKey) return null
    const response = await window.electronAPI.reclusterConditionGraph({ cacheKey, seedNode, options })
    if (!response?.success) {
      set({ error: response?.error || 'Unable to compute clusters for seed node' })
      return response
    }

    set((state) => ({
      reclusterByGraph: {
        ...state.reclusterByGraph,
        [cacheKey]: response,
      },
    }))
    return response
  },

  hydrateReclusterByGraph: (reclusterState = {}) => {
    const safeState = reclusterState && typeof reclusterState === 'object' ? reclusterState : {}
    set({ reclusterByGraph: safeState })
  },

  clearResults: () => set({ results: null, conditionGraphs: {}, reclusterByGraph: {}, error: null }),
}))
