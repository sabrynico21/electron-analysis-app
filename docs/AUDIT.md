# AUDIT — Fedeltà della traduzione Python di SPIEC-EASI (metodo MB)

Data: 2026-10-12
Fase: **1 — Audit (sola lettura)**. Nessun file di codice esistente è stato modificato.

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
`subsample_fraction=0.8`, `stars_threshold=0.05`, `sym_rule='or'`, `random_state=17`,
`require_samples_gt_taxa=False`. Rispetto ai default di R (tutti allineati in Fase 3):
- `n_lambdas=20` ✔ (=R), `lambda_min_ratio=1e-3` ✔ (=R), `n_subsamples=20` ✔ (=`rep.num`);
- `subsample_fraction=0.8` ✔ (=`.ratcheck` per n≤144; per n>144 R userebbe `10√n/n` mentre l'app forza 0.8);
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

### 2.D Warm start vs cold start nel solver MB (scoperta in Fase 2)

Evidenza dai golden (dataset ds4_r, p=10, n=50, senza zeri):

- **Passi deterministici tutti verdi**: CLR `~2e-15`, λ.max `0`, griglia `~1e-16`, coefficienti MB a λ fisso `~1e-13` con **0 edge di differenza**, per **tutti** gli 8 casi.
- **StARS con indici di R iniettati**: identico per `ds1`, `ds2` (diff `≤6e-17`); **diverso** per `ds3`, `ds4` (diff fino a `4e-2`).
- Verifica diretta a λ minimo (ds4_r, subcampione 1): le **adiacenze per-subcampione** di Python coincidono con R se R risolve **un solo λ** (cold) — `0` edge di differenza, β entro `3e-12`. Ma i `merge`/instabilità memorizzati da `spiec.easi` provengono dalla **risoluzione dell'intero path in una sola chiamata huge** (warm), che differisce dalla cold di `max|ΔP|≈0.4`.
- Conferma che R usa il path warm: `trace(pulsar)` mostra che `spiec.easi` passa a pulsar `lambda = path completo`; per ogni subcampione `sparseiCov` risolve tutto il path.

Causa: i dati CLR hanno una dipendenza lineare esatta (le colonne sommano a 0) e, per λ piccoli, l'obiettivo lasso è **piatto/mal condizionato**. Verificato che Python è all'ottimo (valore obiettivo entro `3e-8` da una soluzione sklearn a tolleranza `1e-14`), quindi **non è un bug dell'obiettivo**, ma la superficie ammette soluzioni quasi-ottime con **supporti diversi**; warm e cold ne scelgono una diversa. Conseguenza: l'instabilità StARS può cambiare e, per ds4, il **λ ottimo selezionato cambia** (Python `grid[11]≈1.8e-2` vs R `grid[13]≈8.9e-3`).

## 3. Riepilogo dei verdetti "DIVERSO" (da correggere/decidere in Fase 3)

| Punto | Natura | Impatto |
|-------|--------|---------|
| **warm start vs cold start nel path MB (§1/27, §2.D)** | strategia solver | **λ ottimo e rete finale possono differire** nei casi mal condizionati (ds3, ds4) |
| `nlambda` default 25 vs 20 (§1/5) | default libreria + app | griglia λ più densa di R |
| `lambda.min.ratio` default 1e-2 vs 1e-3 (§1/4) | default libreria + app | path λ più corto; possibile λ ottimo fuori griglia |
| `n_subsamples` default 50 vs 20 (§1/13) | default libreria | più ripetizioni di R (l'app lo corregge a 20) |
| selezione subcampioni / seed (§1/15) | RNG | insiemi di indici diversi ⇒ StARS diverso (iniettabile) |
| bordo `== thresh` nella scelta indice ottimo (§1/18) | logica di confronto | off-by-one solo in uguaglianza esatta |
| floor `tol≥1e-4` (§1/10) | convergenza | impossibile essere più stretti di R |
| `max_iter` ridotto nell'app (§1/10) | convergenza | possibili stime non certificate nel regime n≤10 & p≥200 |
| guardia n>p in libreria (§1/21) | comportamento | errore dove R procede (l'app lo disattiva) |
| path λ non decrescente accettato (§1/8) | validazione input | input che R rifiuterebbe |
| pre-filtro colonne costanti nell'app (§1/22) | comportamento | R fallirebbe invece di scartare |

## 4. Riepilogo dei verdetti "NON IMPLEMENTATO"

- Sostituzione moltiplicativa degli zeri (§1/2) — estensione Python assente in R (default non interessato).
- Spaziatura lineare `lambda.log=FALSE` (§1/7).
- Metodo MB non usa glmnet: il ramo sklearn `Lasso` (§1/24) non ha corrispondenza R ed è fuori dal percorso di default.

## 5. Verifiche numeriche (Fasi 2–3)

Golden R a 17 cifre, 4 dataset × 2 config (`ds1..ds4` × `r`/`app`), indici StARS di R iniettati.

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

## 6. Conclusione (dopo Fase 3)

La pipeline Python MB è ora una traduzione **fedele e verificata** del percorso R: CLR, λ.max, griglia log,
solver huge (warm-started path, ordine/stop/pruning/KKT identici), simmetrizzazione `or`, formula StARS,
cummax, regola di selezione `.starsind`, refit finale. Tutti i passi risultano **EQUIVALENTI/IDENTICI** nei test.

Divergenze rimaste (tutte **intenzionali e documentate**, nessuna altera l'output rispetto a R nella pratica):
1. **Subsampling RNG** (§1/15): RNG diverso; non riproducibile bit-a-bit, mitigato iniettando gli indici di R.
2. **Deviazioni numeriche dichiarate** (§2.A): `np.dot` vs cicli scalari, nessun clamp di S in [−1,1], `hit_max_iter` non esposto.
3. **Scelte intenzionali**: floor `tol≥1e-4`; `max_iter` ridotto a 3000 dall'app in `n≤10 & p≥200` (§1/10);
   ramo sklearn rimosso → `NotImplementedError` (§1/24); guardia `n>p` in libreria (§1/21); pre-filtro colonne costanti nell'app (§1/22).
4. **NON IMPLEMENTATO**: sostituzione moltiplicativa zeri (§1/2); spaziatura lineare `lambda.log=FALSE` (§1/7).

Nessun passo è **NON VERIFICABILE** quanto a sorgente: tutti letti in R/huge/pulsar e verificati numericamente.
