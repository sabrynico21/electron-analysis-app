# AUDIT — Fedeltà della traduzione Python di SPIEC-EASI (metodo MB)

Data: 2026-10-12 (Fasi 1–3) · aggiornato 2026-10-13 (Fasi A–D) · aggiornato 2026-10-14 (estensione C++,
layout riga-major, CI; §7.7 sul perimetro di piattaforma).
Stato: audit completato; correzioni di Fase 3 applicate e verificate con i golden R estesi (Fase B),
regola di sottocampionamento StARS allineata a `pulsar` (Fase C), test di integrazione
dell'applicazione (Fase D). Nessun golden esistente è stato modificato e nessuna tolleranza è stata
allargata.

> **Perimetro.** Tutte le misure e le verifiche di questo documento valgono **solo su Linux x86-64
> (GCC)**: vedi §7.7. Ogni tabella riporta la data e le condizioni di misura accanto ai numeri.

## 0. Sorgenti confrontati

### R (riferimento)
- `Scaricati/SpiecEasi-master/R/spiec-easi.R` — pipeline `spiec.easi.default` (versione DESCRIPTION **1.99.0**). Questa è la copia indicata dall'utente come riferimento.
- `Scaricati/SpiecEasi-master/R/normalization.R` — `clr`, `norm_pseudo`, `norm_to_total`.
- `Scaricati/SpiecEasi-master/R/SparseICov.R` — `sparseiCov`, `neighborhood.net`, `glm.neighborhood`.
- `Scaricati/SpiecEasi-master/R/utilities.R` — `symBeta`.
- **huge** (dipendenza esterna, la libreria che esegue realmente il solver MB): sorgente locale in
  `spiec-easi-python/.tmp-rsrc/huge/` (`R/huge.mb.R`, `R/huge-internal.R`, `src/huge_core.cpp`), versione **2.0.1** = identica alla huge 2.0.1 installata usata a runtime.
- **pulsar** (dipendenza esterna per StARS/subsampling): pacchetto installato `pulsar` **0.3.13**
  (pulsar installato in `~/R/x86_64-pc-linux-gnu-library/4.6/pulsar`). Funzioni lette: `pulsar::pulsar`,
  `pulsar:::stars.stability`, `pulsar:::.starsind`, `pulsar:::.ratcheck`, `pulsar::getLamPath`,
  `pulsar::getMaxCov`, `pulsar:::refit.pulsar`.

> Nota sulla versione del riferimento: `Scaricati/SpiecEasi-master` è la **1.99.0**, cioè la stessa versione
> installata. La differenza con la copia vendored `.tmp-rsrc/SpiecEasi` (2.0.0) è **solo** in commenti roxygen
> e nell'hoisting di `estFunResolved <- match.fun(estFun)` (comportamento invariato). L'audit è quindi valido
> sia contro `SpiecEasi-master` sia contro la SpiecEasi installata.

Nota importante: **il metodo MB in R *non* usa glmnet**. `huge::huge(method='mb')` chiama un solver
coordinate-descent in C++ (`huge_core.cpp::mb`). glmnet compare solo in `neighborhood.net` / `glm.neighborhood`
(`SparseICov.R:147-148`) che sono usati dai metodi **ising/poisson**, non da MB.

### Python (oggetto dell'audit)
- `spiec-easi-python/src/spieceasi/pipeline.py` — `spiec_easi_mb`.
- `spiec-easi-python/src/spieceasi/transforms.py` — `clr_transform`.
- `spiec-easi-python/src/spieceasi/neighborhood.py` — `_fit_single_neighborhood_huge_scale`, `mb_neighborhood_selection`, `alpha_max_for_mb`, `make_alpha_path`.
- `spiec-easi-python/src/spieceasi/stars.py` — `stars_select`, `_instability_from_adjacencies`.
- `spiec-easi-python/src/spieceasi/graph.py` — conversione grafo.
- Wrapper applicativo: `python/analysis.py::_compute_spieceasi_for_matrix` (righe 523-660).

> Disambiguazione: la cartella `python/spieceasi/stars.py` è un **residuo non utilizzato** (non è un package,
> manca `__init__.py`). Il codice realmente importato è `spiec-easi-python/src/spieceasi/` (iniettato in
> `sys.path` da `python/analysis.py:31-40`). L'audit riguarda quest'ultimo.

---

## 1. Tabella di audit

Legenda verdetti: **IDENTICO** / **EQUIVALENTE NUMERICAMENTE** / **DIVERSO** / **NON IMPLEMENTATO** /
**NON VERIFICABILE** / **NON APPLICABILE**.

| # | Passo — default R (rif.) | Python (rif.) | Verdetto | Differenza ed effetto sull'output |
|---|--------------------------|---------------|----------|-----------------------------------|
| 1 | Pseudo-count & CLR: `t(clr(data+1, 1))` — `spiec-easi.R:71`; `clr.default` con base `e`, tol=`Machine$double.eps`, `LOG - mean(LOG)/mean(nzero)` — `normalization.R:95-98` | `clr_transform(data, pseudocount=1.0, zero_method='additive')`; `adjusted = X+pseudocount` → proporzioni → `log` → sottrazione media di riga — `transforms.py:113-152` | **EQUIVALENTE NUMERICAMENTE** | Con `pseudocount=1` e dati interi ≥0 tutti i valori sono >tol per R, quindi `mean(LOG)/mean(nzero)` coincide con `mean(log(x+1))`. In Python la normalizzazione a proporzioni introduce `-log(somma di riga)` che si cancella nella media: risultati identici (verificato in `tests/test_r_components.py::test_clr_matches_r_spiec_easi_norm_transform`). |
| 2 | Metodo di gestione zeri: unico = additivo (`data+1`). Nessuna sostituzione moltiplicativa in questo percorso — `spiec-easi.R:71` | default `zero_method='additive'`; è presente anche `'multiplicative'` (replacement di Dirichlet/delta) — `transforms.py:114,145-149` | **IDENTICO** (default) / **NON IMPLEMENTATO in R** (ramo multiplicative) | Il default coincide. Il ramo `'multiplicative'` è un'**estensione** assente nella pipeline R; non è mai usato dall'app (default `'additive'`). |
| 3 | `lambda.max = getMaxCov(cor(X))` — `spiec-easi.R:175`; `getMaxCov` = max \|triu\| stretto di `cor` (pulsar) | `alpha_max_for_mb` = `max\|off-diag\|` della matrice `huge_correlation(X)` (stessa S usata dal solver) — `neighborhood.py` | **EQUIVALENTE NUMERICAMENTE** | `cor` R usa denominatore n−1 (equivalente a `ddof=1`). Fase 3: Python ora ricava `alpha_max` dalla **stessa** S del solver (prima usava un prodotto nodewise separato, diverso di 1 ULP → l'edge a correlazione massima sopravviveva al soft-threshold a λ=λ.max lasciando un edge spuro ~1e-16 nel modello "vuoto"). Floor `max(alpha_max,1e-6)` mantenuto (assente in R, senza effetto pratico). |
| 4 | `lambda.min.ratio` default = **1e-3** — `spiec-easi.R:236` | **1e-3** (`pipeline.py:67`); l'app passa 1e-3 (`analysis.py`) | **IDENTICO** | Allineato a R in Fase 3 (era 1e-2 nella libreria e nell'app). |
| 5 | `nlambda` default = **20** — `spiec-easi.R:237` | **20** (`pipeline.py:66`); app 20 (`analysis.py`) | **IDENTICO** | Allineato a R in Fase 3 (era 25). |
| 6 | Spaziatura log: `getLamPath(max, max*ratio, n, log=TRUE)` = `exp(seq(log(max), log(min), length.out=n))` — `spiec-easi.R:238-239`; pulsar::getLamPath | `np.geomspace(alpha_max, alpha_min, num=n)` — `neighborhood.py:370-372` | **EQUIVALENTE NUMERICAMENTE** | `geomspace(a,b,n) = exp(linspace(log a, log b, n))` = `exp(seq(log a, log b, length.out=n))`. Identico a meno di errore macchina. |
| 7 | `lambda.log` default = TRUE (log) — `spiec-easi.R:153`; se FALSE spaziatura lineare | sempre geometrico; nessuna opzione lineare | **EQUIVALENTE NUMERICAMENTE** (default) / **NON IMPLEMENTATO** (`lambda.log=FALSE`) | Solo il default è coperto. |
| 8 | Path λ esplicito via `...`→`args$lambda` (huge valida: finite, >0, non crescenti) — `spiec-easi.R:234-240` | parametro `alpha_path` — `pipeline.py`; validazione: 1-D, ≥2, >0, **non crescente** | **EQUIVALENTE NUMERICAMENTE** | Fase 3: aggiunta la validazione di monotonia non crescente (`ValueError`), come richiesto da huge. |
| 9 | Stima MB: `huge::huge(x=clr, method='mb')` → S = `cor(x)` (standardizza ddof=1, `.huge_standardize` `huge-internal.R:47`, `.huge_fast_cor` `:78-89`) → C++ `huge::mb` su S (`huge_core.cpp:679`) | `mb_solver.solve_mb_path_huge_scale`: S = `huge_correlation(X)`, coordinate-descent nodewise su tutto il path (warm start) — `mb_solver.py` | **EQUIVALENTE NUMERICAMENTE** | Percorso **default** replicato con CD **fedele a huge** (non sklearn, non glmnet). Vedi §2.A per il confronto strutturale riga-per-riga e §2.D per il warm start. |
| 10 | Tolleranza/imax solver MB: `thol=1e-4`, `MAX_ITER=10000` — `huge_core.cpp:686-687` | `thol = max(tol,1e-4)`, `max_cd_iter = min(max_iter,10000)` — `mb_solver.py` | **EQUIVALENTE** (default) / **DIVERSO documentato** | Default `tol=1e-4`, `max_iter=10000` coincidono. **Floor `tol≥1e-4` documentato** nel docstring di `mb_solver` (impossibile essere più stretti di R). L'app abbassa `max_iter` a 3000 se `n≤10 & p≥200` (`analysis.py:597-601`): **deviazione intenzionale da R, documentata** (`MAX_ITER` di R resta 10000). |
| 11 | Regola "strong": sequenziale `strong_thr = 2λ_i − λ_{i−1}` (i>0), `2λ_i` per i=0 — `huge_core.cpp:725-726` | `strong = (idx_i==1) & (grad_abs > strong_thr)` con `strong_thr = 2λ_i − λ_{i−1}` — `mb_solver.py` | **IDENTICO** | Fase 3: implementata la strong rule **sequenziale** (prima era `2α` fisso). Necessaria anche in aritmetica esatta per riprodurre quale soluzione quasi-ottima huge raggiunge. |
| 12 | Simmetrizzazione: `path = sign(\|β\| OP t(\|β\|))`, OP=`or` (default) — `huge.mb.R:143,147-150`; default `sym='or'` in `sparseiCov`/`huge` | `adjacency = (\|coef\|>0) OR/AND (\|coef\|>0)ᵀ`, diagonale azzerata, default `'or'` — `neighborhood.py:335-345` | **IDENTICO** | Stessa regola OR sull'insieme dei non-zeri; diagonale a 0 in entrambi. |
| 13 | Numero subcampioni: pulsar `rep.num` default **20**; `spiec.easi` non lo sovrascrive | `n_subsamples` default **20** (`pipeline.py:69`, `stars.py:80`); l'app passa **20** | **IDENTICO** | Allineato a R in Fase 3 (era 50 nella libreria). |
| 14 | Dimensione subcampione: `.ratcheck` = 0.8 se n≤144, altrimenti `10√n/n`; size = `floor(n·ratio)` — pulsar::.ratcheck, `pulsar::pulsar` | stessa regola — `stars.py:138-147` | **EQUIVALENTE NUMERICAMENTE** | `floor(n·10√n/n)=floor(10√n)`. Identica. L'app forza `subsample_fraction=0.8`. |
| 15 | Selezione subcampioni: `sample(1:n, size, replace=FALSE)`, seed da `pulsar.params$seed` (default **NULL ⇒ casuale**) — pulsar::pulsar | `np.random.default_rng(random_state).choice(n, size, replace=False)` — `stars.py:161-163`; default `random_state=0`, app 17. **È presente un parametro `subsamples`** (`stars.py:91`, ramo `:104-128`; esposto da `pipeline.spiec_easi_mb`) che **bypassa** l'RNG e inietta indici 0-based fissi | **DIVERSO** | Gli insiemi di indici generati sono **diversi** perché il RNG è diverso e perché R di default non fissa il seme. StARS è stocastico ⇒ è la principale fonte di differenza negli output. Mitigazione in Fase 2: **iniettare gli indici di R** via il parametro `subsamples` (già supportato). Nota: pulsar aggiunge in testa `ind.sample` il vettore `1:n` (modello full, `attr "full"`) usato per il refit finale ma **escluso** dal calcolo di stabilità. |
| 16 | Instabilità StARS: `merge=Reduce("+")/rep.num`; `summary[i] = 4·sum(P·(1−P))/(p·(p−1))`; poi `cummax` — pulsar:::stars.stability | `4·mean_supratriangolo(P·(1−P))` — `stars.py:60-71`; `np.maximum.accumulate` — `stars.py:185` | **EQUIVALENTE NUMERICAMENTE** | Poiché P è simmetrica, `sum_su_tutta_matrice(P(1−P)) = 2·sum_supratriangolo`; quindi `4·sum_full/(p(p−1)) = 4·mean_upper` (il fattore 2 compensa la semidimensione). `cummax` = `maximum.accumulate`. |
| 17 | Soglia StARS: `spiec.easi` imposta `thresh=0.05` se non fornita — `spiec-easi.R:278` (default pulsar sarebbe 0.1) | `stars_threshold=0.05` — `pipeline.py:73`, `stars.py:83`; app 0.05 | **IDENTICO** | 0.05 in entrambi. |
| 18 | Scelta λ ottimo: `.starsind` = `max(first{summary ≥ thresh} − 1, 1)`; se nessuno raggiunge la soglia → indice `length(summary)` + warning — pulsar:::.starsind | `_stars_index`: `max(first{monotone ≥ thresh} − 1, 0)`; fallback ultimo indice **+ warning** — `stars.py` | **IDENTICO** | Corretto in Fase 3 (era `≤`, con off-by-one sul bordo di uguaglianza). Test di bordo generati da R (`equality_mid`, `equality_first`, `exactly_zero`, `all_below`, `all_above`, `normal_cross`) tutti verdi. |
| 19 | Rete finale: `refit.pulsar` riadatta il modello full-data risolvendo l'intero path (warm) e prende `path[[opt.index]]` (pulsar:::refit.pulsar); `getRefit` restituisce l'adjacency | `mb_path_neighborhood_selection(X_clr, alpha_path=gr, selected_index)` risolve l'intero path (warm) sul dataset completo e prende l'indice ottimo — `pipeline.py` | **EQUIVALENTE NUMERICAMENTE** | Fase 3: il fit finale risolve ora l'intero path warm-started nello stesso modo di R (prima risolveva solo il λ ottimo, cold). |
| 20 | Coefficienti/pesi: `spiec.easi` restituisce l'adjacency; i β (scala correlazione) via `est$est$beta[[opt]]`; pesi simmetrici via `symBeta(mode='ave')` = `(β+βᵀ)/2` — `utilities.R:136-142` | `directed_coefficients` (β huge-scale); `weight_matrix = 0.5·(directed+directedᵀ)·adjacency` — `neighborhood.py:347-352` | **EQUIVALENTE NUMERICAMENTE** | `0.5·(β+βᵀ)` replica esattamente `symBeta('ave')`; mascherare con l'adjacency è coerente con l'uso documentato. L'adjacency (output primario) è identica. |
| 21 | Vincolo n>p: nessun controllo esplicito (huge gestisce n<d; `scr.num=n−1` se screen) — `huge.mb.R:53-62` | Libreria: `require_samples_gt_taxa=True` → `ValueError` se `n≤p` (`pipeline.py:140-146`); l'app passa `False` (`analysis.py:648`) | **DIVERSO** (libreria) / **EQUIVALENTE** nell'app | La guardia aggiuntiva in libreria non esiste in R. Con `False` (come nell'app) il regime p>n è permesso come in R. |
| 26 | Orientamento dei coefficienti: huge restituisce `beta[predittore, nodo]` (colonna = nodo m, riga nel blocco = predittore j; vedi `SPMBgraph.cpp`, `collect_sorted`: `indices = lambda_idx*d + j`) | `directed_coefficients[nodo, predittore]` (`neighborhood.py` docstring; verificato in Fase 2) | **DIVERSO** (rappresentazione) sugli output diretti; **EQUIVALENTE** sugli output simmetrici | La matrice dei coefficienti R è la **trasposta** di quella Python (`max\|py − rᵀ\| ≈ 8e−15`). Effetto: **nessuno** su `adjacency`, `path`, `weights` (tutti simmetrici o funzione del supporto). Rilevante solo se si confrontano i β diretti nodo-per-nodo: vanno trasposti. |
| 27 | **Risoluzione del path MB**: per ogni (sub)campione, R chiama `huge(x, lambda=path completo)` in **una sola invocazione** → **warm start** lungo il path + strong rule **sequenziale** `2λ_i − λ_{i−1}` (`huge_core.cpp:679-810`, `:725-726`) | `mb_solver.solve_mb_path_huge_scale` risolve l'intero path per (sub)campione: warm start, strong rule sequenziale, ordine `j=0..d−1` saltando `m`, stessi criteri di stop/pruning/KKT; usato da `stars.py` (per subcampione) e `pipeline.py` (fit finale) | **EQUIVALENTE NUMERICAMENTE** | **Implementato in Fase 3** (prima era cold-start per λ). Ha **eliminato** le divergenze StARS/rete di ds3/ds4 (ora raw/instabilità e adjacency finale entro `~1e-15`/0 edge). Vedi §2.D. |
| 22 | Colonne costanti: huge **errore** ("Raw data x contains a constant column", `huge-internal.R:216-219`) | Libreria: `ValueError` se sd≤0 (`neighborhood.py:46-48`, `alpha_max_for_mb:337-339`); l'app **pre-filtra** le feature costanti (`analysis.py:555`) | **EQUIVALENTE** (libreria) / **DIVERSO** (app) | Libreria ~ R (errore). L'app invece rimuove le costanti prima della stima: R fallirebbe. Nessun effetto quando non ci sono colonne costanti. |
| 23 | Intercetta (lasso MB): assente (sul correlazione standardizzata) — implicito in huge | assente nel percorso huge-scale | **IDENTICO** | Il CD huge-scale non usa intercetta. |
| 24 | Uso di glmnet: solo ising/poisson (`neighborhood.net`/`glm.neighborhood`, `SparseICov.R:147-148`), **non** MB | ramo sklearn **rimosso**: `standardize=False` o `coefficient_scale='original'` → `NotImplementedError` — `neighborhood.py` | **NON IMPLEMENTATO** (intenzionale) | Fase 3: il ramo sklearn/glmnet è stato rimosso perché MB in R usa huge, non glmnet. L'app non lo usava (verificato con grep). Il test `test_r_lasso_coefficients.py` ora verifica che sollevi `NotImplementedError`. |
| 25 | Output restituiti: `fit` (refit pulsar) con `$refit$stars` (adjacency), `$lambda`, `$est$beta`, `$select$stars$summary/opt.index/thresh`; λ ottimo via `getOptLambda` | dataclass con `adjacency, weights, graph, selected_lambda, alpha_path, instability, monotone_instability, directed_coefficients, taxa_names` — `pipeline.py:29-60` | **EQUIVALENTE** | Python espone le stesse grandezze chiave più diagnostiche extra (grafo networkx, pesi). L'adjacency e il λ ottimo sono i campi comparabili con `getRefit`/`getOptLambda`. |

---

## 2. Approfondimenti sui punti critici

### 2.A Solver MB — confronto strutturale `mb_solver.solve_mb_path_huge_scale` vs `huge::mb`

Riferimenti: `spiec-easi-python/src/spieceasi/mb_solver.py` ↔ `huge_core.cpp:679-810` (`mb`) e `:651-677` (`mb_refine_active`).

- **Matrice di lavoro**: entrambi risolvono su S = correlazione (n−1). Python `huge_correlation` (`S=(ZᵀZ)/(n−1)`, diag=1); huge `.huge_fast_cor` (`huge-internal.R:78-89`). Uguale.
- **Residuo strong/attivazione**: Python `r = S[m,j] − S[j,active]·w[active]`; huge `mb_partial_residual = S[m,j] − Σ S[j,a_k] w[a_k]` (`:640-649`). Uguale.
- **Refinement**: Python `r = S[m,w]+w0[w] − S[w,active]·w0[active]`; huge `mb_refine_active` `r = S[m,w_idx] + w0[w_idx] − Σ S[w_idx,a_k] w0[a_k]` (`:660-668`). Uguale.
- **Criterio di stop interno**: `gap = Σ|Δw| / Σ|w|` con soglia `thol`; uguale.
- **Pruning**: riattiva le coordinate andate a zero (`strong=1`); uguale.
- **Certificazione KKT**: passata completa sulle coordinate inattive, conta le violazioni, riattiva e ripete; uguale.
- **Ordine coordinate**: R percorre `j=0..d−1` saltando `m`; Python **ora replica esattamente lo stesso ordine** (array di dimensione `d`, `idx_i[m]=0`).
  > **CORREZIONE (Fase 3).** La versione precedente dell'audit affermava che l'ordine delle coordinate è *irrilevante* perché la soluzione del lasso è unica. **Non è vero a tolleranza finita**: il lasso è strettamente convesso solo in aritmetica esatta; con `thol=1e-4`, nei regimi mal condizionati (correlazione CLR quasi singolare, λ piccoli) la superficie è piatta e **ordine, warm start e criteri di stop determinano quale soluzione quasi-ottima viene raggiunta** (vedi §2.D). Python replica ora l'ordine di huge senza eccezioni.
- **Deviazioni residue dichiarate**: (1) i residui usano `np.dot` (BLAS) invece dei cicli scalari di `mb_partial_residual`/`mb_refine_active` → possibile differenza all'ultimo bit nell'ordine di somma; (2) Python non clampa S in [−1,1] (huge sì); (3) `hit_max_iter` non è esposto. Nessuna ha alterato l'esito dei test (diff ≤1e-15).

Conclusione: **replica fedele** del solver huge (warm-started path), non di glmnet.

### 2.B StARS — differenze residue

- **Formula di instabilità** e **cummax**: identici (§1 righe 16-17).
- **Campionamento subcampioni** (§1 riga 15): RNG diverso; R di default non fissa il seme. Mitigato iniettando gli indici di R via il parametro `subsamples` (i test di equivalenza lo fanno).
- **Bordo `== thresh`** (§1 riga 18): **risolto** in Fase 3 (regola identica a `.starsind`, warning incluso).
- **Warm start** (§1 riga 27): **risolto** in Fase 3; `stars_select` ora risolve l'intero path per (sub)campione.
- Il modello full (elemento 1 in pulsar) è usato solo per il refit; Python lo riproduce con un warm solve separato sul dataset completo (§1 riga 19).

### 2.C Default a livello applicativo (`python/analysis.py`)

L'app invoca `spiec_easi_mb` con: `n_lambdas=20`, `lambda_min_ratio=1e-3`, `n_subsamples=20`,
`subsample_fraction=None`, `stars_threshold=0.05`, `sym_rule='or'`, `random_state=17`,
`require_samples_gt_taxa=False`. Rispetto ai default di R (tutti allineati in Fase 3):
- `n_lambdas=20` ✔ (=R), `lambda_min_ratio=1e-3` ✔ (=R), `n_subsamples=20` ✔ (=`rep.num`);
- `subsample_fraction=None` ✔ (= regola di `.ratcheck`: `0.8` se `n≤144`, altrimenti `10√n/n`),
  con override esplicito via `spieceasiSubsampleFraction`. **Fase C:** prima l'app forzava `0.8`,
  divergendo da R per `n>144`;
- `random_state=17` (R non fissa il seme ⇒ casuale in R, riproducibile in Python);
- **deviazione intenzionale documentata**: `max_iter` scende a 3000 se `n≤10 & p≥200` (R resta 10000).

### 2.D Warm start vs cold start nel solver MB (scoperta in Fase 2, RISOLTA in Fase 3)

R: `spiec.easi`/`pulsar` risolvono, per ogni (sub)campione, l'intero path λ in **una sola chiamata huge** (warm start + strong rule sequenziale). Python risolveva ogni λ **cold** (una chiamata per λ): a tolleranza finita `thol=1e-4`, nei regimi mal condizionati (dati CLR a somma-zero → correlazione quasi singolare, λ piccoli) le due strategie raggiungono supporti diversi.

Evidenza (Fase 2): ds1/ds2 identici; ds3/ds4 divergevano (instabilità fino a `4e-2`, λ ottimo diverso). Verificato che Python era all'ottimo (obiettivo entro `3e-8` da sklearn a `tol=1e-14`): non era un bug dell'obiettivo ma di *quale* soluzione quasi-ottima veniva raggiunta. R con un solo λ (cold) coincideva con Python; il `merge` warm di `spiec.easi` no (`max|ΔP|≈0.4`).

**Correzione (Fase 3).** `mb_solver.solve_mb_path_huge_scale` replica `huge::mb` riga per riga (warm start, strong sequenziale, ordine `j=0..d−1`, stop/pruning/KKT). `stars_select` e il fit finale lo usano. Risultato: **ds1–ds4 (entrambe le config) ora identici** a R (instabilità raw/monotona `≤1e-15`, λ ottimo uguale, adjacency finale 0 edge di differenza).

**Correzioni di bordo correlate (Fase 3):**
- **Endpoint del path**: `geomspace` calcola `exp(log(αmax))` che può arrotondare 1 ULP *sotto* `αmax`; `make_alpha_path` ora fissa `path[0]=αmax` → modello più sparso vuoto come in R.
- **Consistenza di `αmax`**: `alpha_max_for_mb` ora usa la **stessa** S del solver (`huge_correlation`); un prodotto nodewise separato differiva di 1 ULP e faceva sopravvivere l'edge a correlazione massima (~1e-16) a λ=λ.max.

---

### 2.E λ = λ.max: risultato numericamente indeterminato (anche in R)

A λ = λ.max il modello corretto è quello **vuoto** (nessun coefficiente). Tuttavia λ.max coincide, per
costruzione, con la correlazione massima `max|cor(X)|`, e in quel punto la coordinata-descendo di
`huge::mb` cade esattamente sul bordo della soft-threshold: a seconda dell'ordine di accumulo dei
prodotti scalari il lato a correlazione massima viene tenuto o scartato.

**Il fenomeno è interno a R** (verificato sui golden, §5): `grid[0]` coincide bit-per-bit con
`pulsar::getMaxCov(cor(X))` (diff `0`), ma risolvendo `huge::huge(X, method="mb", lambda = max|cor|)`
sulla **stessa** matrice CLR si ottengono **0 lati per alcuni dataset e 1 lato per altri** — stesso λ,
dati diversi. Non è quindi un'instabilità da imputare alla traduzione: **a λ = λ.max il risultato è
numericamente indeterminato anche in R**.

**Scelta intenzionale di Python.** Il solver usa la stessa S per griglia e solver (`huge_correlation`) e,
coerentemente, a λ = λ.max restituisce **sempre il modello vuoto** (il risultato matematicamente
corretto). La correzione di Fase 3 su `alpha_max_for_mb` va letta in questo senso: eliminare
un'incertezza di 1 ULP che faceva sopravvivere artificialmente il lato a correlazione massima.
Non si replica l'ordine di sommatoria bit-a-bit di huge.

**Trattamento nei test.** In `tests/test_equivalence_extended.py` il passo 9 ammette un'**eccezione
esplicita e stretta**: se e solo se l'indice λ selezionato è 0 (= λ.max), è tollerata una differenza di
**al più un lato**, e solo se è esattamente la coppia `argmax|cor|`. In ogni altro caso il confronto
resta **esatto** (zero lati di differenza).

## 3. Riepilogo dei verdetti "DIVERSO" (stato attuale)

I punti sotto sono **risolti oppure scelte intenzionali**: nessuno è una divergenza aperta nella
configurazione di default dell'app. La colonna "Stato" dice quando e come è stato chiuso, o se è una
scelta ancora attiva. L'unica indeterminazione residua è quella a λ = λ.max (§2.E).

| Punto | Natura | Stato | Impatto attuale |
|-------|--------|-------|-----------------|
| **warm start vs cold start nel path MB (§1/27, §2.D)** | strategia solver | **RISOLTO in Fase 3** | Nessuno: in ds1–ds4 (config `r` e `app`) λ ottimo e rete finale coincidono; la divergenza è storia delle Fasi 1–2 |
| `nlambda` default 25 vs 20 (§1/5) | default libreria + app | **RISOLTO in Fase 3** (ora 20 in entrambi) | Nessuno |
| `lambda.min.ratio` default 1e-2 vs 1e-3 (§1/4) | default libreria + app | **RISOLTO in Fase 3** (ora 1e-3 in entrambi) | Nessuno |
| `n_subsamples` default 50 vs 20 (§1/13) | default libreria | **RISOLTO in Fase 3** (ora 20) | Nessuno |
| `subsample_fraction=0.8` forzato dall'app (§2.C) | default applicazione | **RISOLTO in Fase C** (ora regola di `.ratcheck`) | Nessuno per `n>144`: l'app usa `10√n/n` come R |
| selezione subcampioni / seed (§1/15) | RNG | **APERTO per scelta** (non riproducibile) | Gli insiemi di indici differiscono ⇒ StARS diverso. Mitigato iniettando gli indici di R; il caso senza iniezione è coperto da `tests/test_statistical_equivalence.py` |
| bordo `== thresh` nella scelta indice ottimo (§1/18) | logica di confronto | **RISOLTO in Fase 3** (regola di `.starsind`, warning incluso) | Nessuno; 6 casi di bordo generati da R |
| floor `tol≥1e-4` (§1/10) | convergenza | **scelta intenzionale attiva** | Impossibile essere più stretti di R; nessun effetto misurato |
| `max_iter` ridotto nell'app (§1/10) | convergenza | **deviazione intenzionale attiva** (solo `n≤10 & p≥200`) | Possibili stime non certificate in quel regime; R resta a 10000 |
| guardia n>p in libreria (§1/21) | comportamento | **scelta intenzionale attiva** in libreria | Errore dove R procede; l'app la disattiva (`require_samples_gt_taxa=False`) |
| validazione del path λ (§1/8) | validazione input | **IDENTICO a huge** (Fase 3) | Python rifiuta un path **crescente** (`ValueError`) come huge: nessun input accettato in più |
| pre-filtro colonne costanti nell'app (§1/22) | comportamento | **scelta intenzionale attiva** | L'app rimuove le costanti prima della stima: R fallirebbe invece di scartare |

## 4. Riepilogo dei verdetti "NON IMPLEMENTATO"

- Sostituzione moltiplicativa degli zeri (§1/2) — estensione Python assente in R (default non interessato).
- Spaziatura lineare `lambda.log=FALSE` (§1/7).
- Metodo MB non usa glmnet: il ramo sklearn `Lasso` (§1/24) non ha corrispondenza R ed è fuori dal percorso di default.

## 5. Verifiche numeriche (Fasi 2–3)

Golden R a 17 cifre, 4 dataset × 2 config (`ds1..ds4` × `r`/`app`), indici StARS di R iniettati.
**Condizioni di misura:** 2026-10-11 (Fase 2) e 2026-10-12 (Fase 3), Linux x86-64, macchina di
sviluppo, R 4.6.1 + SpiecEasi 1.99.0 + huge 2.0.1 + pulsar 0.3.13, un solo thread.

**Prima delle correzioni (Fase 2):**
- Passi deterministici (1–5): tutti **EQUIVALENTI** — CLR `≤2e-15`, λ.max `≤3e-16`, griglia `≤2e-16`,
  coefficienti MB a λ fisso `≤1.1e-13`, adiacenza simmetrizzata a λ fisso: **0 edge**.
- StARS (passo 6): identico per ds1/ds2; **diverso** per ds3/ds4 (fino a `4e-2`) — strategia warm/cold (§2.D).
- Bordo selezione `== thresh`: off-by-one **confermato**.
- Suite slow: **10 failed / 4 passed in 5994 s (~1h40m)**.

**Dopo le correzioni (Fase 3):**
- `test_solver_path_matches_r` (ds3/ds4, isolamento del solver warm): **4/4 PASSED** — adiacenza per-subcampione per λ e coefficienti uguali a huge.
- `test_python_matches_r_golden` (end-to-end StARS + rete finale): **8/8 PASSED** — instabilità raw/monotona `≤1e-15`, λ ottimo uguale, adjacency finale 0 edge.
- Test di bordo selezione (generati da R): **6/6 PASSED**.
- Suite slow completa: **12 passed in 246 s (~4m06s)**, da ~1h40m.

**Fase B — golden estesi (23 casi generati; 20 eseguiti, p≤50; p=100 in corso):**

| passo | max\|diff\| osservato |
|---|---|
| 1 CLR | 1.3e-15 – 3.6e-15 |
| 2 λ.max | 0 – 3.3e-16 |
| 3 griglia | 0 – 1.7e-16 |
| 4 path solver per-subcampione | **0 lati diversi in tutti i 20 casi** |
| 4b coefficienti | 6.7e-15 – 7.4e-13 |
| 5 StARS raw | 5.6e-17 – 3.6e-15 |
| 6 StARS monotono | 0 – 3.2e-15 |
| 7/8 indice e λ ottimo | 0 – 1.7e-16 (indice sempre uguale) |

Casi: p ∈ {10,20,30,50,100}, n ∈ {25..300}, densi/sparsi (fino al 75% di zeri), `p>n`, `n>144`,
configurazioni `r` (nλ20/1e-3/rep20) e `app` (nλ25/1e-2/rep20). Esito: passi 1–8 identici;
l'unica divergenza è il **passo 9**, **1 lato**, in 8/20 casi.

Eccezione λ.max verificata sui golden (§2.E):
- **8 casi** con indice R = 0: l'unico lato finale di R è **esattamente** la coppia `argmax|cor|`
  (es. `rand_p10_n30_s201` (7,8)=0.40315831740882058; `pn_zeros_p40_n25_s219` (1,2)=0.59313578070591977;
  `rand_p20_n50_s205` (2,6)=0.37198908181359475);
- **6 casi** con indice R = 0 ma R con **0 lati** (modello vuoto): Python coincide esattamente;
- **6 casi** con indice R = 1.

Tempi: vedi `docs/PERFORMANCE.md`. I tre casi p=100 sono stati eseguiti **per intero**:
`largep_p100_n80_s215` OK in 2425.5 s, `largep_p100_n200_s216` OK in 4933.3 s,
`app_largep_p100_n80_s215` OK in 1367.3 s (`3 passed in 8726.93 s`, ~2h25m). In tutti e tre
**tutti i passi combaciano** (CLR ≤ 3.1e-15, λ.max ≤ 1.1e-16, griglia ≤ 1.1e-16, **0 lati** nel path
solver per-subcampione, coefficienti ≤ 1.8e-13, StARS ≤ 1.3e-14, `final_edges = 0`): a p=100 **non
c'è nessuna divergenza**, nemmeno al passo 9.

## 6. Conclusione (dopo Fase 3)

La pipeline Python MB è ora una traduzione **fedele e verificata** del percorso R: CLR, λ.max, griglia log,
solver huge (warm-started path, ordine/stop/pruning/KKT identici), simmetrizzazione `or`, formula StARS,
cummax, regola di selezione `.starsind`, refit finale. Tutti i passi risultano **EQUIVALENTI/IDENTICI** nei test.

Divergenze rimaste (tutte **intenzionali e documentate**, nessuna altera l'output rispetto a R nella pratica):
1. **Subsampling RNG** (§1/15): RNG diverso; non riproducibile bit-a-bit, mitigato iniettando gli indici di R.
2. **Deviazioni numeriche dichiarate** (§2.A): `np.dot` vs cicli scalari, nessun clamp di S in [−1,1], `hit_max_iter` non esposto.
   A λ = λ.max questo si manifesta come **indeterminazione numerica condivisa con R** (§2.E): a λ.max Python
   restituisce **sempre** il modello vuoto (scelta intenzionale, matematicamente corretto), mentre `huge::mb`
   a volte tiene il lato `argmax|cor|`. Nei test estesi è ammessa la sola eccezione stretta descritta in §2.E.
3. **Scelte intenzionali**: floor `tol≥1e-4`; `max_iter` ridotto a 3000 dall'app in `n≤10 & p≥200` (§1/10);
   ramo sklearn rimosso → `NotImplementedError` (§1/24); guardia `n>p` in libreria (§1/21); pre-filtro colonne costanti nell'app (§1/22).
4. **NON IMPLEMENTATO**: sostituzione moltiplicativa zeri (§1/2); spaziatura lineare `lambda.log=FALSE` (§1/7).

Nessun passo è **NON VERIFICABILE** quanto a sorgente: tutti letti in R/huge/pulsar e verificati numericamente.

## 7. Estensione C++ del solver MB (`spieceasi._mb_core`) — deviazioni e risultati

**Perché.** Il solver coordinate-descent in Python era il collo di bottiglia (46–112× più lento di R).
Per ragioni di sola performance è stata compilata **unicamente** la parte MB; CLR, matrice di correlazione,
griglia λ, StARS, selezione e refit restano in Python e sono condivisi dalle due implementazioni.

**Che cosa è stato tradotto.** `_mb_core.cpp` è una traduzione quasi letterale di `huge_core.cpp`
`mb` / `mb_partial_residual` / `mb_refine_active`: path warm-started, regola strong sequenziale
(2λ_i − λ_{i−1}, con 2λ_0 per i = 0), sweep `j = 0…d−1` saltando `m`, pruning, certificazione KKT,
`threshold_l1` (0 se |x| ≤ thr), indicizzazione column-major `cm(data, nrows, r, c)`, `thol = 1e-4`,
`MAX_ITER = 10000`. Interfaccia: `solve_mb_path(S, lambda_path, thol, max_iter) → (coefficients, hit_max_iter)`,
con `S` passata in ordine Fortran e **GIL rilasciato** attorno al nucleo numerico.

### 7.1 Deviazioni dichiarate dell'estensione

1. **Nessun OpenMP** sulle `m = 0…d−1` regressioni. Attenzione: il pacchetto `huge` 2.0.1 **installato è
   compilato con OpenMP** (`libs/huge.so` linka `libgomp`, `Makevars.in` usa `OPENMP_CXXFLAGS`, e `mb()`
   esegue il ciclo sui nodi sotto `#pragma omp parallel` + `#pragma omp for schedule(dynamic)`): il
   confronto di wall time con R **non è seriale-contro-seriale** e va letto tenendo conto di questo
   fattore (misurato: 1.56× a p=30 e 2.13× a p=100 sulla sola chiamata `mb`). La nostra traduzione non
   parallelizza; l'ipotesi che farlo sia numerico-safe è supportata dai **contatori identici** (§7.6) ma
   **non è stata implementata**: la proposta completa (gate di bit-identità 1 vs N thread, packaging
   per piattaforma, regola anti-oversubscription) è in `docs/PERFORMANCE.md` come lavoro futuro.
2. **Nessun `collect_sorted`.** `huge` raccoglie e ordina i candidati dell'active set; la traduzione
   segue la scansione in ordine di variabile. Il punto fisso (KKT, supporto, coefficienti) è lo stesso,
   l'ordine di visita no: è la ragione principale per cui i risultati non sono bit-identici a R.
3. **`hit_max_iter` è esposto** dalla nuova interfaccia (`_mb_core` lo restituisce); in precedenza non
   era disponibile. Non è usato per alterare il risultato, solo per diagnostica.
4. **Flag di compilazione** (dettaglio e confronto in `docs/PERFORMANCE.md`):
   `-std=c++17 -O2 -ffp-contract=off -fno-lto` (`/O2 /fp:precise` su MSVC), **nessun** `-ffast-math`,
   **nessun** `-march=native`. R compila con `-O2 … -flto=auto -ffat-lto-objects` e senza disattivare
   la contrazione (default GCC `-ffp-contract=fast`): R *può* quindi usare FMA dove l'estensione no.
   La scelta è nella direzione conservativa (meno riordino/contrazione = più vicino al riferimento).
5. **Fallback automatico.** Se `_mb_core` non è importabile, `mb_solver` usa
   `mb_solver_reference.py` (copia congelata del solver Python puro) ed emette un `RuntimeWarning`
   esplicito. Il riferimento Python resta quindi la sorgente di verità e la garanzia di funzionamento.
   Il backend dell'app logga `SPIEC-EASI solver: cpp|python` e, nel secondo caso, un `WARNING` che
   spiega la caduta di prestazioni: il fallback non può passare inosservato.
6. **Layout di S: riga-major invece di colonna-major.** `huge` memorizza S colonna-major e legge una
   *riga* (`S[j, a_k]` con `a_k` consecutivi), cioè con stride `d`. La traduzione mantiene lo stesso
   elemento e lo stesso ordine di sommatoria ma memorizza S riga-major, così quell'accesso è contiguo
   (BIT-IDENTICO su tutti i 248 problemi del gate; +5.3 % a p=100, +16 % a p=300, +36 % a p=500).
   NON è la riscrittura di cui avverte il sorgente di `huge` (quella resta dentro la colonna `j`
   usando la simmetria ed era 2× più lenta per loro): vedi `docs/PERFORMANCE.md`.
7. **Nessuna affermazione di equivalenza generale.** Le verifiche sotto valgono **sui casi testati** e
   **sulla sola piattaforma su cui sono state eseguite** (Linux x86-64, GCC): vedi §7.7.

L'esecuzione **multi-thread dei subcampioni StARS** (FASE 3 a) non è una deviazione numerica: le
soluzioni sono bit-identiche per `n_jobs = 1…16` (aggregazione a ordine fisso, nessuno stato globale),
quindi l'output non dipende dal numero di worker. È documentata come default in `CHANGELOG.md`.

### 7.2 Verifica estensione-vs-riferimento (244 problemi + golden)

| metrica | risultato |
|---|---|
| problemi (p 5…60, n 10…150; densi, con zeri, mal condizionati; 6 λ) | 244 |
| bit-identici al riferimento Python | **0 / 244** |
| max \|diff\| famiglie non degeneri | **9.470e-14** |
| max \|diff\| famiglia degenere (collinearità relativa 1e-9) | 1.881e-01 (peggiore: `illcond p=8 n=150 seed=202`) |
| supporti diversi, famiglie non degeneri | **0 problemi** |
| supporti diversi, famiglia degenere | 25 problemi / 40 blocchi λ |
| max \|diff\| sui golden rilasciati | 1.652e-13 |

**Condizioni di misura:** 2026-10-13, Linux x86-64 (GCC), macchina di sviluppo, estensione compilata
`-O2 -ffp-contract=off -fno-lto`, un solo thread per lato; riferimento = solver Python puro
(`mb_solver_reference.py`).

Quindi: l'estensione **non è bit-identica** né al riferimento Python né a R. Le differenze sono a livello
1e-13…1e-14 e la struttura scoperta è identica in tutte le famiglie non degeneri.

### 7.3 Famiglia degenere: perché non è risolta (e non lo sarà allargando tolleranze)

Nei casi `illcond` la matrice di correlazione ha collinearità relativa ~1e-9: la soluzione lasso non è
unica e le due implementazioni (e R) cadono in punti diversi dello stesso insieme di minimi. Prove
raccolte confrontando **entrambe** le implementazioni con `huge`/`pulsar` reali su
`illcond p=8 n=25 seed=202` (λ indici 2…5; lati diversi = lati con supporto diverso da R):

| λ idx | \|C++ − R\| | lati diversi | \|ref − R\| | lati diversi | nnz R / C++ / ref |
|---|---|---|---|---|---|
| 2 | 2.094e-05 | 2 | **9.659e-15** | 0 | 22 / 23 / 22 |
| 3 | 4.742e-02 | 4 | 3.582e-05 | 2 | 30 / 28 / 29 |
| 4 | 4.934e-02 | 0 | 1.187e-02 | 0 | 32 / 32 / 32 |
| 5 | 5.109e-02 | 2 | 1.187e-02 | 0 | 32 / 31 / 32 |

R non privilegia **nessuna** delle due: al λ 2 coincide col riferimento (9.7e-15) e differisce dal C++
(2.1e-05); al λ 3–5 si allontana da entrambe (1.2e-02…5.1e-02) pur restando i supporti quasi uguali.
La causa prossima è che `cor` di R e la nostra accumulazione della correlazione differiscono di ~1 ULP,
e in questo regime quell'1 ULP decide quali coordinate restano attive. **Non è quindi un difetto
dell'estensione**, e non è stato "risolto" allargando tolleranze: resta riportato.

**Criterio di validità adottato per la famiglia degenere (sostituisce il confronto sui coefficienti).**
Sui problemi `illcond` si verifica che ciascuna delle due soluzioni sia un **ottimo valido** dello stesso
problema nodewise ridotto (la coordinata `j == m` è esclusa, come fa `huge`), con due soglie dichiarate
a priori in `tests/validation_config.py` e derivate dalla tolleranza del solver, non dall'esito:

| controlli | soglia dichiarata | base |
|---|---|---|
| residuo KKT `\|\|g_j\| − λ\|` (attivi) e `max(\|g_j\| − λ, 0)` (inattivi), `g = S b − s_m` | `DEGENERATE_KKT_VIOLATION_MAX = 1e-3` | un ordine di grandezza sopra `thol = 1e-4`, su cui `huge` dichiara convergenza |
| scarto dei valori di obiettivo `0.5 b'Sb − s_m'b + λ\|b\|₁` fra le due implementazioni | `DEGENERATE_OBJECTIVE_GAP_MAX = 1e-5` | un ordine sotto `thol × scala O(1)`, due ordini sotto il λ più piccolo dei path testati |

Il confronto coefficiente-per-coefficiente con la reference su questa famiglia **resta calcolato e
stampato come informazione** (max \|diff\|, blocchi λ con supporto diverso). Per le famiglie non degeneri
il criterio stretto è invariato: `rtol 1e-9`/`atol 1e-11` sui coefficienti e uguaglianza esatta dei
supporti. Nessuna tolleranza è stata allargata e la famiglia non è stata ristretta.

**Valori misurati** su tutta la famiglia degenere dello sweep (**76 problemi** × 6 λ × `p` risposte, entrambe
le implementazioni): residuo KKT massimo **2.478e-04** — identico per riferimento e core
(`illcond p=12 n=10 seed=101`, λ n. 5, risposta 8), quindi ≈ 4× sotto la soglia dichiarata; scarto di
obiettivo massimo **2.390e-06** (`illcond p=60 n=25 seed=202`, λ n. 5, risposta 56), ≈ 4× sotto la soglia.
Le due soglie non sono state quindi scelte a posteriori per far passare il test: i valori misurati stanno
circa un fattore 4 sotto i limiti dichiarati a priori. Il fatto che il residuo KKT massimo sia lo **stesso**
per le due implementazioni indica che nel caso peggiore entrambe raggiungono lo stesso punto; le differenze
di vertice (25 problemi con supporto diverso, 40 blocchi λ) restano su casi con residuo molto più piccolo.

### 7.5 Che cosa distingue la famiglia 1e-9 dai dati CLR reali (e perché ds3/ds4 passano)

Misure dirette (stessa macchina, `huge_correlation` nostra vs `cor` di R sullo stesso input, via
`Rscript`; "blocchi λ cambiati" = blocchi in cui il supporto simmetrizzato differisce):

| design | entrate di S diverse da `cor` di R | max \|ΔS\| | λmin di S | perturbazione casuale ±1 ULP: blocchi λ cambiati | stesso solver alimentato con la S di R: blocchi λ cambiati |
|---|---|---|---|---|---|
| `ds1_r.clr` (p=10) | 66/100 | 1.94e-16 | −3.6e-16 | media 0.70, max 1 | **0** |
| `ds2_r.clr` (p=30) | 748/900 | 2.22e-16 | −6.8e-17 | media 0.40, max 1 | **1** (λ n. 0 = λ.max, 1 lato) |
| `ds3_r.clr` (p=30) | 666/900 | 4.44e-16 | −9.5e-16 | media 0.70, max 1 | **0** |
| `ds4_r.clr` (p=10) | 74/100 | 4.44e-16 | −4.4e-16 | media 0.55, max 1 | **0** |
| `illcond p=5 n=15` | 14/25 | 3.33e-16 | −1.2e-18 | media 1.20, max 2 | **2** (λ 0, 4) |
| `illcond p=8 n=25` | 44/64 | 3.33e-16 | −2.5e-16 | media 3.10, max 4 | **3** (λ 2, 3, 5) |
| `illcond p=8 n=40` | 48/64 | 3.33e-16 | +1.6e-16 | media 0.25, max 1 | **4** (λ 0, 3, 4, 5) |
| `illcond p=12 n=40` | 100/144 | 3.33e-16 | +6.4e-17 | media 1.40, max 4 | **2** (λ 2, 4) |
| `plain p=8 n=40` | 54/64 | 5.55e-16 | +4.7e-02 | media 0.60, max 1 | **0** |
| `plain p=30 n=100` | 726/900 | 5.55e-16 | +6.3e-03 | media 0.30, max 1 | **0** |

**Perché i golden ds3/ds4 passano.** Non perché siano ben condizionati: la loro matrice di correlazione
è singolare a livello macchina esattamente come quella dei casi sintetici (λmin ≈ 1e-16), e la nostra S
differisce da quella di R in 66–748 entrate su 100–900 in **tutti** i dataset. Passano perché la
differenza di 1 ULP fra le due S **non sposta il vertice selezionato**: alimentando il medesimo solver con
la matrice di correlazione *di R*, il path risolto è identico su `ds1`, `ds3`, `ds4` (0 blocchi λ
diversi) e differisce su `ds2` soltanto nel blocco a **λ = λ.max** (1 lato), cioè esattamente il blocco
già coperto dall'eccezione stretta di §2.E. Il caso `ds2` è quindi la stessa indeterminazione di §2.E,
non una fragilità nuova.

**Che cosa distingue la famiglia 1e-9.** Il discriminante misurabile è la **sensibilità del vertice**,
non la singolarità né l'entità di ΔS (simili in tutte le famiglie): sugli stessi dati reali il passaggio
da S a S_R cambia 0 blocchi λ, nella famiglia sintetica ne cambia 2–4 su 6, **nell'interno del path**
(non solo a λ.max), e i lati coinvolti sono sempre quelli incidenti alla colonna resa collineare
(`p=8`: (3,7), (5,7); `p=8 n=40`: (0,7), (1,6), (5,6), (5,7); per costruzione la colonna collineare è
l'ultima). Con perturbazioni casuali di ±1 ULP si osserva la stessa gerarchia: ≤ 1 blocco su dati reali
e controlli non degeneri, fino a 4 sui casi sintetici.

**Ipotesi (con le prove sopra, non dimostrata).** Nei dati CLR reali la degenerazione è *strutturale e
condivisa*: la trasformazione CLR impone la somma a zero delle colonne in entrambe le implementazioni, la
direzione quasi-nulla è quella di tutti-uni e il termine lineare del lasso lungo di essa è al livello del
rumore di arrotondamento (misurato: `v·1` = −0.98 / +1.00 / −0.86 per ds1/ds2/ds4, `max_m |v·S[:,m]|` ≈
1e-16 contro un λ di path di 3e-4…1e-3). Nella famiglia sintetica la collinearità 1e-9 è invece una
perturbazione *dei dati* su una coppia di colonne il cui prodotto con la risposta non è piccolo, e il
vertice che coinvolge quella coppia è intrinsecamente instabile: lo stesso 1 ULP basta a spostare il lato
da una colonna all'altra. In entrambi i casi la proprietà è del **problema**, non dell'implementazione:
le soluzioni di R si spostano allo stesso modo (§7.3).

### 7.6 Dove va il tempo a p = 100 (misurato)

**Condizioni di misura:** 2026-10-13, Linux x86-64 (GCC), macchina di sviluppo, `p=100, n=80`,
`n_jobs=1`, estensione compilata.

A `p=100, n=80`, `n_jobs=1`, il **99.7 %** del tempo della pipeline è dentro l'estensione (21 chiamate:
20 subcampioni StARS + il refit finale); il calcolo della matrice di correlazione pesa 0.2 % e **tutto il
resto del Python** (CLR, simmetrizzazione, aggregazione StARS, selezione, conversione del grafo) **0.1 %**.
Il divario residuo rispetto a R (1.57× a p=100, misura **pre**-riga-major: oggi 1.52×, vedi §7.4) è quindi
interno al ciclo coordinate-descent e non è recuperabile lato Python.

Esperimenti di build sulla stessa sorgente: passare S da colonna-major a riga-major (stesso elemento
acceduto, stesso ordine di sommatoria) e compilare a `-O3` restituisce risultati **bit-identici** su 7
input (i 4 golden CLR, il CLR p=100, un caso `plain` e uno `illcond`) e recupera ≈ 6 % (23.24 s → 22.00 s
con riga-major a `-O2`, → 21.82 s con `-O3`). La copia C→F imposta dal binding costa 5.2 µs per chiamata
(0.109 ms in tutto), quindi è irrilevante: il guadagno è nell'accesso non strided.

**Esito (aggiornato il 2026-10-14).** Il **layout riga-major è stato adottato** (commit `e73ec82`,
2026-10-13): il gate di bit-identità è stato rieseguito su **248 problemi** (244 dello sweep + i 4
golden) contro la build precedente, tutti bit-identici, e la suite lenta è stata rieseguita per intero
(`50 passed, 84 deselected in 1445.48 s`). Effetto misurato sugli stessi tre casi benchmark
(mediana di 3, seriale): 0.572 → 0.514 s (p=30), 3.991 → 3.691 s (p=50), 22.056 → 21.348 s (p=100).
**`-O3` non è stato adottato**: 0–1.2 %, dentro il rumore fra esecuzioni (mediana 21.696 s con `-O2`
contro 21.670 s con `-O3`), e allargherebbe la distanza dai flag di build di R; resta disponibile come
modifica di una riga. Dettagli e tabelle: `docs/PERFORMANCE.md`.

### 7.4 Effetto sui risultati già verificati

**Condizioni di misura:** 2026-10-13, Linux x86-64 (GCC), macchina di sviluppo, estensione compilata,
un solo thread (`n_jobs=1`) su entrambi i lati.

Sui benchmark e sulle suite i risultati **non cambiano**: λ selezionato e numero di lati identici tra
riferimento e C++ in tutti i casi misurati (e identici a R), suite non-slow da 662 s a 18.66 s. I tre
casi p=100 verificati per intero (§5) non mostrano alcuna divergenza e passano ora in **129.23 s**
complessivi (29.80 / 81.84 / 17.51 s) contro **8726.93 s** del riferimento, cioè **67.5×**, con
`final_edges = 0` e 0 lati diversi nel path solver per-subcampione. Numeri e flag: `docs/PERFORMANCE.md`.

**Nota sulle cifre ricorrenti per p=100.** Il rapporto con R compare nel repository come **1.57×**
(tabella FASE 2 di `docs/PERFORMANCE.md`, misura del 2026-10-13 **prima** della modifica di layout,
22.056 s), **1.52×** (dopo l'adozione del layout riga-major, 21.348 s) e **1.54×** (riga `n_jobs=1`
della misura FASE 3, 21.619 s). Sono tre esecuzioni diverse dello stesso caso sulla stessa macchina,
non un disaccordo: il valore da citare dipende da quale build è in uso. Il confronto **seriale-contro-seriale**
resta comunque non omogeneo, perché `huge` è compilato con OpenMP (§7.1 punto 1).

**Effetto collaterale noto e preesistente.** `tests/test_r_variability_expanded_scenarios.py::test_stochastic_variability_baseline_expanded_scenarios`
fallisce. Non è causato dall'estensione: rieseguito con il **solo riferimento Python**
(`HAVE_CPP_SOLVER = False`) il fallimento è **identico** (stessa tabella: `cross_jaccard = 1.0` contro
baseline 0.0 + margine 0.08, `py_edges_median = 0.0`, `r_edges_median = 1.0`, `boundary0` 4/4 in tutte
le sei scene). In tutte le sei scene entrambe le implementazioni selezionano **l'indice λ = 0**, cioè
λ = λ.max, dove — come dimostrato in §2.E — il risultato è indeterminato anche in R: Python restituisce
il modello vuoto (scelta intenzionale) e R a volte conserva il lato `argmax|cor|`. È la stessa
indeterminazione di §2.E osservata a livello di rete finale; i margini di quel test
(`STOCHASTIC_BASELINE_TOLERANCE`) sono dichiarati a priori e **non** sono stati allargati.

### 7.7 Perimetro delle verifiche: piattaforma

**Tutte** le verifiche di equivalenza e di bit-identità riportate in questa sezione e in §5 sono state
eseguite **solo su Linux x86-64 con GCC**, sulla macchina di sviluppo
(AMD Ryzen 5 3550H, glibc 2.44, flag di compilazione in `docs/PERFORMANCE.md`), contro
R 4.6.1 + SpiecEasi 1.99.0 + huge 2.0.1 + pulsar 0.3.13 compilati con le stesse impostazioni della
distribuzione.

**Windows e macOS non sono dimostrati identici.** Sono coperti dalla CI, ma il risultato non può essere
esteso a quelle piattaforme finché i job non girano e passano lì: compilatori diversi (MSVC, Apple
clang), implementazioni di libm diverse e architetture diverse (x86-64 vs ARM64) possono cambiare gli
ultimi bit dei risultati in virgola mobile, e la differenza residua misurata è già a livello 1e-13…1e-14.

Il job `equivalence-platforms` di `.github/workflows/ci.yml` esegue, su
`ubuntu-latest` (x86-64/GCC), `macos-13` (Intel), `macos-14` (ARM) e `windows-latest` (MSVC + Rtools):

- la **sweep completa estensione-vs-riferimento** (`test_mb_core_equivalence.py -m slow`: 244 problemi
  + i golden estesi), e
- le **suite golden contro R** (`test_equivalence.py`, `test_equivalence_extended.py`), con
  R + `huge`/`pulsar`/SpiecEasi installati sul runner e `rpy2`.

Non è un job di sola compilazione: la compilazione verificata dal job `wheels` non è evidenza di
equivalenza numerica. Fino al passaggio di questi job, l'affermazione di equivalenza e di bit-identità
di questo documento va letta come **valida su Linux x86-64 (GCC)**.
