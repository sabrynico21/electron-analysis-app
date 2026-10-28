# Public microbiome dataset candidates for equivalence validation

**Purpose:** identify candidate independent data sources before generating
additional reference fixtures. This is an inventory, not approval to redistribute
the source data. Verify the original data terms and attribution before committing
any data or derived fixtures.

## Existing baseline

| Dataset | Source | Current use |
|---|---|---|
| QMP | `SPRING-master/data/QMP.rda` | Existing public dataset; fixtures cover 106×20 and 20×30 subsets |

## Candidates inspected in the local R environment

| Candidate | Public source / provenance | Availability and observed profile | Status |
|---|---|---|---|
| American Gut Project Round 1 (`amgut1.filt`) | `SpiecEasi` R package data; help attributes source to the American Gut Project / Human Food Project | Available locally through `SpiecEasi` 1.99.0 (GPL-3); 289 samples × 127 OTUs; 30.64% zero entries; no all-zero samples | Second-study source used in the fixtures. Provenance recorded in `DATASET_PROVENANCE.md`; redistribution **not cleared** (`redistribution_cleared: false`) |
| GlobalPatterns / Caporaso et al. 2011 (`phyloseq::GlobalPatterns`) | Caporaso et al., PNAS 2011, DOI `10.1073/pnas.1000080107`, PMCID `PMC3063599`; article says data deposited as SRA Study `SRA012609.1` | Available locally through phyloseq 1.56.0 (AGPL-3); 26 samples × 19,216 OTUs as packaged, 79.07% zero entries, 228 all-zero taxa, no all-zero samples. Core 26×15 and 26×40 fixtures use the most prevalent varying taxa and have 0% zeros. A prevalence-spread 26×40 variant is retained as a CCLasso reference-failure stress input. | Independent-study source used in the fixtures. Provenance recorded in `DATASET_PROVENANCE.md`; redistribution **not cleared** (`redistribution_cleared: false`). Fixture matrices omit metadata/tree |
| American Gut Project Round 2 (`amgut2.filt.phy`) | Same `SpiecEasi` data documentation / American Gut Project | 296 samples × 138 taxa; 34.60% zeros; 2 all-zero samples; no sample-ID overlap with the Round 1 matrix | Useful robustness case, but same study/project as Round 1, so do **not** count as an independent study in the paper's dataset count. Exclude or explicitly handle the two all-zero samples using an R-matched preprocessing rule |
| Human Microbiome Project 16S (`hmp216S`, listed as `hmp216S (hmp2)`) | Listed by `SpiecEasi` as Human Microbiome Project 2 data | Listed in the installed `SpiecEasi` dataset index but supplied by the separate `hmp2` package, which is not installed in this environment | Strong independent-study candidate. Obtain the versioned `hmp2` source/dependency, inspect dimensions, zero/library-size profile, metadata, license, and preprocessing before selecting fixture cases |

American Gut Round 1 and GlobalPatterns are the two selected independent study
sources beyond QMP. The two American Gut rounds have disjoint sample IDs locally,
but their shared study design means they are not treated as independent datasets
for the primary breadth claim. HMP 16S remains a reserve candidate; its `hmp2`
data package is not installed locally.

The fixtures now use deterministic cases from both selected sources and retain
sample/taxon IDs, source dimensions, zero fraction, library-size range, and source
label. Provenance and redistribution status are recorded per fixture and audited
in `DATASET_PROVENANCE.md`. **Disposition (2026-10-20): neither American Gut Round
1 nor GlobalPatterns is cleared for redistribution.** The derived matrices are
retained only as internal CI fixtures (`redistribution_cleared: false`) and must
not be published or shipped until written permission is obtained or the cases
are replaced with synthetic data.

## Required checks before accepting a candidate

1. Cite the source publication and repository/package version; archive the data
   access URL and retrieval date.
2. Confirm data-use/redistribution terms, especially for human-derived data.
   Do not commit raw participant-level data without explicit permission.
3. Record source dimensions, taxa/sample identifiers, zero fraction, all-zero
   rows/columns, library-size quantiles, and metadata exclusions.
4. Write deterministic preprocessing code shared by the R reference generator
   and Python fixture consumer. Do not hand-select samples/taxa after observing
   method outputs.
5. Select representative `p < n` and `p > n` cases and sparse/less-sparse
   profiles under the acceptance contract. Store only permitted derived fixtures,
   source hashes, preprocessing scripts, and provenance manifests.
