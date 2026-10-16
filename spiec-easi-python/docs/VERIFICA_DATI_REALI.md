# VERIFICA R↔Python sui dati reali e su dati sintetici — disegno pre-registrato

Questo documento è **pre-registrato**: definizioni dei casi, metriche e soglie di
lettura sono fissate **prima** di eseguire la batteria. Le soglie qui dichiarate
**non vengono modificate** dopo aver visto i risultati. I risultati vengono
riportati in appendice (sezioni A e successive) senza toccare le regole.

Obiettivo: una batteria di verifica R↔Python **sistematica e realistica** sui
dati reali (`genere_table.txt` + `sample_metadatab.txt`) e su dati sintetici,
**solo a parità di subcampioni** (gli indici di `pulsar` sono esportati da R e
iniettati in Python). Nessuna modifica a golden, tolleranze, eccezioni o codice
di produzione: **solo script, test slow e documentazione**.

---

## 0. Contesto e regole generali

- Stack di riferimento: **SpiecEasi 1.99.0 / huge 2.0.1 / pulsar 0.3.13**
  (R 4.6.1, come verificato con `packageVersion`). Python: numpy 2.5.3,
  scipy 1.18.1, pandas 3.0.6 nel venv `spiec-easi-python/.venv`.
- Configurazione di default per **tutti** i casi (salvo i sintetici dove
  indicato): `nlambda=20`, `lambda.min.ratio=1e-3`, `rep.num=20`,
  soglia StARS `0.05`, regola di simmetrizzazione `"or"`, `max_iter=10000`,
  frazione di subcampione da `pulsar::.ratcheck`
  (`0.8` se `n<=144`, altrimenti `floor(10*sqrt(n))/n`, dimensione
  `floor(ratio*n)`).
- Dati reali: matrice conteggi `1135 genera × 96 campioni` (`genere_table.txt`,
  trasposta: le righe sono i generi). Metadati: 96 campioni, `Condition` con
  **48 `Treated`** e **48 `Not_Treated`**.
- **Subcampioni condivisi**: gli indici (0-based) usati da R sono salvati e
  iniettati in Python via il parametro `subsamples`. R non è riproducibile con
  il RNG di numpy, quindi l'unico protocollo valido per il confronto è
  iniettare gli stessi indici. Ogni caso ha **10 seed di R** (10 estrazioni di
  indici di subcampione), dichiarati.
- **Indici di λ**: R è **1-based** (`opt.index` parte da 1), Python è
  **0-based**. La conversione (`python_index = r_index - 1`) è esplicita e
  coperta da un **test unitario** dedicato.
- **Primo passo divergente**: se un caso diverge, si riporta il primo passo
  divergente con contesto (subcampione, indice λ, differenza). **Non si
  allargano le tolleranze.**
- Il caso con `p=1135` **non** viene usato per confronti di tempo finché non è
  chiarito il punto F.

### Metriche di lettura (soglie a priori, non modificabili)

| Passo | Grandezza | Soglia |
| --- | --- | --- |
| 1 | CLR | `rtol=1e-12, atol=1e-12` |
| 2 | `λ.max` | `rtol=1e-12, atol=1e-15` |
| 3 | griglia λ | `rtol=1e-12, atol=1e-15` |
| 4a | supporto solver per subcampione/λ | 0 lati diversi (eccezione λ.max: ≤1 lato = coppia `argmax|cor|`, cfr. `docs/AUDIT.md` §2.E) |
| 4b | coefficienti solver | `rtol=1e-6, atol=1e-8` |
| 5 | curva StARS grezza | `rtol=1e-6, atol=1e-8` |
| 6 | curva StARS cumulata | `rtol=1e-6, atol=1e-8` |
| 7 | indice λ scelto | uguaglianza esatta (dopo conversione 1-based/0-based) |
| 8 | valore λ scelto | `rtol=1e-9, atol=0` |
| 9 | rete finale | 0 lati diversi (eccezione λ.max come 4a) |
| validità | KKT (massima violazione) | `≤ 1e-3` |
| validità | scarto di obiettivo vs miglior solver | `≤ 1e-5` |

Per la **classificazione** (sezione D) si usa inoltre: curva StARS
`max|diff| ≤ 1e-10` per l'esito "identico".

### Non vacuità dei casi

Per ogni caso si riporta **quanti seed su 10 scelgono un λ interno**
(`opt_index > 0`). Se quasi tutti scelgono `λ.max` (indice 0), il caso è
**dichiarato degenere** e non contribuisce alle conclusioni sulla selezione.

---

## A. Prevalenza e dimensioni (stima prima di eseguire)

Sui 96 campioni (`all`) e sui sottogruppi `Treated` / `Not_Treated` (48 ciascuno)
si calcola, per ciascuna soglia di prevalenza, il numero di generi `p` (generi
presenti in almeno `soglia×100%` dei campioni del gruppo), la frazione di zeri
della sotto-matrice risultante e `n`. La **definizione di prevalenza** di un
genere è la frazione di campioni del gruppo con conteggio `> 0`.

### Regola di scelta delle soglie (dichiarata a priori)

Dalla tabella A si sceglie una soglia il cui `p` cada nel bersaglio:

- **R1, R2** → `p ≈ 50–150`;
- **R3** → `p ≈ 150–250`.

Regola operativa: **R1 = soglia più bassa con `p ≤ 150`** (caso realistico più
grande), **R2 = soglia che porta `p` nella metà inferiore di [50,150]**, **R3 =
soglia più bassa con `150 ≤ p ≤ 250`**. Se le soglie nominali di B
(≥50%, ≥20%, ≥10%) danno `p` fuori bersaglio, si usano le soglie ottenute con
questo criterio e la discrepanza viene **dichiarata nel report** (non si riducono
i parametri in silenzio).

### Stima dei tempi in R

I tempi di R sono stimati per ogni caso da una **micro-calibrazione misurata**
(un fit `mb` su path completo, pochi subcampioni) e dalla scalatura attesa
`≈ O(p² · n · nlambda · rep.num)` con il numero di lati effettivi. Se una stima
supera **30–60 minuti in R**, il caso si **ferma** e si propongono alternative
(senza ridurre i parametri in silenzio).

---

## B. Casi della batteria

Ogni caso realistico/limite è eseguito con **10 seed di R** (10 estrazioni di
subcampione); gli indici di λ e di subcampione sono esportati da R e iniettati
in Python.

### Casi sui dati reali (n = 96 campioni)

| ID | Gruppo | Soglia prevalenza | n | Bersaglio `p` |
| --- | --- | --- | --- | --- |
| R1 | tutti | da A (regola) | 96 | 50–150 |
| R2 | tutti | da A (regola) | 96 | 50–150 |
| R3 | tutti | da A (regola) | 96 | 150–250 |
| R4 | Treated | ≥20% nel sottogruppo | 48 | (da A) |
| R5 | Not_Treated | ≥20% nel sottogruppo | 48 | (da A) |
| R6 | tutti | ≥20% | 96 | 415 |

Nota: R4/R5 usano la **stessa regola** di selezione applicata al sottogruppo.

**R6** è il punto intermedio tra denso e sparso richiesto esplicitamente:
96 campioni, soglia ≥20% → `p=415`, `zero_fraction ≈ 0.38`. Prima di eseguire i
10 seed si **misura** il tempo reale di R su **un** seed (la stima di A.3 in quel
range è estrapolata e instabile): se supera ~30–60 minuti per seed, il caso si
ferma e si propongono alternative **senza ridurre i parametri in silenzio**.

### Nota densità / sparsità (aggiunta dopo l'approvazione di A)

I casi **R1–R5 sono dataset densi** (frazione di zeri ≈ 3–9%) e quindi
sostengono **solo** l'affermazione "equivalenza su dati densi". La **sparsità
reale** del dataset (~74% di zeri, cfr. L2) è coperta **solo** da:

- **L2** (100 generi casuali, zero_frac ≈ 0.66–0.79),
- **L3** (≥5% di prevalenza, `p=691`, zero_frac ≈ 0.585),
- **R6** (≥20%, `p=415`, zero_frac ≈ 0.38) — ponte intermedio,
- i **casi sintetici** S1–S4 (frazione di zeri dichiarata nel modello).

Nessuna conclusione va estesa oltre questi casi.

### Casi limite

| ID | Definizione |
| --- | --- |
| L1 | 24 campioni casuali su 96, soglia ≥20%; **10 estrazioni diverse** con seed dichiarati |
| L2 | 100 generi casuali senza filtro; **10 estrazioni diverse**, inclusa quella del vecchio caso D |
| L3 | soglia ≥5%, `p` limitato dai tempi (da A) |
| L4 | 100 generi più abbondanti (per abbondanza media totale), controllo |

### Casi sintetici (verità nota)

Modello generativo dichiarato: conteggi multinomiali (depth noto) da un modello
di Abbondanza Log-Normale `mu = exp(scale(Z·Wᵀ))`, dove `Z` è latente
`n×k`, `W` è `p×k`, `k = max(1, p//4)`, e una matrice di **precisione sparsa
nota** (banda/grafo a lati fissi) su cui si applica CLR. Segnale calibrato per
ottenere almeno un λ interno scelto.

| ID | p | n | Segnale | Note |
| --- | --- | --- | --- | --- |
| S1 | 50 | 100 | debole | |
| S2 | 50 | 100 | medio | |
| S3 | 50 | 100 | forte | |
| S4 | 60 | 200 | medio | regola `n>144` (`ratcheck`) |

Per S1–S4 si dichiarano matrice di precisione, `k`, `depth`, frazione zeri e il
seed; la verità (supporto) è nota e si confronta la rete finale con essa a
titolo descrittivo (concordanza con la verità), **non** come criterio di parità
R↔Python.

### Casi patologici (confronto di comportamento)

| ID | Definizione |
| --- | --- |
| P1 | due generi identici o proporzionali |
| P2 | un genere presente in un solo campione |
| P3 | un campione con pochissime letture |

Per P1–P3 si riporta se **R e Python danno lo stesso risultato o lo stesso
errore** (messaggio di errore comparabile), non necessariamente gli stessi
valori.

---

## B.1 Decisioni a priori sull'esecuzione (fissate prima di eseguire)

### Ordine e blocchi

1. **(a)** `R1–R5, L1, L2, L4, S1–S4, P1–P3` con **10 seed** ciascuno (seed
   `1001…1010`);
2. **(b)** `R6` con **5 seed**; **(c)** `L3` con **5 seed**;
3. **(d)** punto F (`p=1135`, senza confronti di tempo).

Dopo **(a)** ci si ferma e si riporta il **report intermedio** (classificazione D
dei casi eseguiti). L'esecuzione avviene **in background**, a blocchi, salvando i
CSV grezzi dopo ogni (caso, seed); lo script è **idempotente** (salta le coppie
(caso, seed) già presenti nel CSV) e riprendibile.

### Estensione a 10 seed di R6/L3 (regola a priori)

R6 e L3 partono con **5 seed**. Si estendono a **10 seed** *solo se l'esito dei
primi 5 è **ambiguo***, dove "ambiguo" significa **almeno una** delle:

- classificazione (sezione D) **non uniforme** tra i 5 seed (mescolanza di
  `identico`/`equivalente`/`diverso`);
- l'**indice λ scelto** (`opt.index`) **cambia** tra i 5 seed;
- comparsa di `r_error`/`py_error`/divergenza in una parte dei seed ma non in
  tutti.

Se tutti e 5 i seed danno lo **stesso esito**, ci si ferma a 5 e si riporta
l'**intervallo di Wilson 95% su n=5**, **dichiarando esplicitamente che il
campione è piccolo** e che non si estende.

### Nota sui tempi (Python vs R)

I tempi Python sono misurati con `n_jobs=-1` (**thread**, parallelismo sui
subcampioni); `pulsar` in R esegue invece **seriale**. Il rapporto tra i tempi
**non misura l'efficienza del solver**: confronta due modalità di esecuzione
diverse (Python parallelo vs R seriale). Va letto solo come ordine di grandezza
del costo della batteria.

---

## C. Misure per ogni caso e seed

1. **Passi deterministici**: CLR, `λ.max`, griglia λ, dimensione del subcampione.
2. **Solver isolato sull'intero path** (tutti i 20 λ) sul dataset completo:
   supporto e coefficienti R vs Python.
3. **Solver per subcampione** (i 20 subcampioni di R): lati diversi per λ e
   coefficienti, con attenzione ai **λ profondi** (grandi indici).
4. **Curva StARS** (grezza e cumulata), λ scelto, rete finale (lati, lati
   diversi), pesi.
5. **Dove il solver diverge**:
   - (a) criterio di **validità già dichiarato** (KKT `1e-3`, scarto obiettivo
     `1e-5`) applicato a **entrambe** le soluzioni;
   - (b) diagnostica **"S di R"**: si rifà il calcolo Python alimentando il
     solver con la **matrice di correlazione di R** e si riporta se la
     divergenza scompare.
6. **Non vacuità**: per ogni caso, quanti seed scelgono un λ interno
   (indice `> 0`).

---

## D. Classificazione decisa a priori

Per ogni caso, esito tra:

- **"identico entro tolleranza"**: stessa rete (0 lati diversi), stesso λ
  scelto, curva StARS `max|diff| ≤ 1e-10`, e tutti i passi deterministici entro
  le tolleranze della tabella §0.
- **"equivalente per validità"**: le soluzioni del solver **differiscono** ma
  **entrambe valide** per KKT (`≤1e-3`) e scarto di obiettivo (`≤1e-5`), **mentre
  rete finale e λ scelto sono uguali**.
- **"diverso"**: rete finale o λ scelto diversi. In tal caso si riporta il
  **primo passo divergente** (subcampione, indice λ, differenza) e l'analisi.

Precedenza: se rete/λ differiscono → "diverso"; altrimenti se tutto entro
`1e-10` → "identico"; altrimenti (solver diverso ma entrambe valido) →
"equivalente per validità".

---

## E. Output

- Tabella per caso con le **frequenze** ("λ uguale in X/10 seed", "rete identica
  in X/10", "solver identico sul path in X/10"), con **intervalli di Wilson**
  (95%) dove ha senso, e una **tabella riassuntiva**.
- **Script riproducibili**, **CSV grezzi** con **SHA-256**, e **test slow** che
  rilegge i risultati e verifica le soglie dichiarate (i test **non** rilanciano
  R).
- Conclusioni **non più forti dei dati**: nessun "equivalente sui tuoi dati" in
  generale, solo **frequenze osservate per quei casi**. Si elencano **dati non
  verificabili** e **incongruenze**.

---

## F. Caso `p=1135` (chiarimento, senza confronti di tempo)

Si chiama in R direttamente `huge::huge(method="mb")` con la **stessa `S`** e gli
**stessi λ** su **un solo subcampione**, sia con `OMP_NUM_THREADS=1` sia con
tutti i thread. Si contano **sweep** e **aggiornamenti di coordinata** in R e in
Python (gli stessi contatori già usati per `p=100`) e si verifica che i λ siano
quelli **reali** del path. Si riporta se il **lavoro è identico**. In caso
negativo, il confronto di tempo precedente viene **scartato**.

---

## Appendice A — Risultati (compilata dopo l'esecuzione della stima)

> Le regole e le soglie della sezione A **non** sono state modificate. Questa
> appendice aggiunge solo i numeri osservati e le conseguenti scelte.

### A.1 Tabella prevalenza / dimensioni

Script: `spiec-easi-python/tests/real_data/estimate_dimensions.py`.
Input (in `spiec-easi-python/tests/real_data/data/`):

- `genere_table.txt` — SHA-256 `4632669c214ad03a927c5cb22adebc0f2e3eb8dd825a2251d45ffcb64810253b`
- `sample_metadatab.txt` — SHA-256 `09d339e3a19bb94c414a869c857de33bf253949fba535d58d9ed4a71b86fe914`

Output grezzo: `spiec-easi-python/tests/real_data/dimensions.json`
(SHA-256 `08e71bbea6727e1a568a91622fc92956cd9a483aeb7bd2fea1b1e89e9ca0fdc8`).

Matrice totale: **1135 generi × 96 campioni**; gruppi: `Treated` 48,
`Not_Treated` 48.

`p` = generi con conteggio `> 0` in almeno `soglia` frazione dei campioni del
gruppo. `zero_fraction` = frazione di zeri nella sotto-matrice selezionata.

**Gruppo `all` (n = 96):**

| soglia | p | zero_fraction |
| --- | --- | --- |
| 0.05 | 691 | 0.5851 |
| 0.10 | 555 | 0.5005 |
| 0.20 | 415 | 0.3796 |
| 0.30 | 342 | 0.2998 |
| 0.40 | 286 | 0.2319 |
| 0.50 | 245 | 0.1785 |
| 0.60 | 207 | 0.1265 |
| 0.70 | 174 | 0.0857 |
| 0.75 | 163 | 0.0729 |
| 0.80 | 150 | 0.0590 |
| 0.85 | 134 | 0.0445 |
| 0.90 | 114 | 0.0297 |
| 0.95 | 82 | 0.0130 |

**Gruppo `Treated` (n = 48):**

| soglia | p | zero_fraction |
| --- | --- | --- |
| 0.05 | 672 | 0.5694 |
| 0.10 | 558 | 0.4960 |
| 0.20 | 432 | 0.3885 |
| 0.30 | 347 | 0.2979 |
| 0.40 | 291 | 0.2295 |
| 0.50 | 249 | 0.1739 |
| 0.60 | 211 | 0.1226 |
| 0.70 | 181 | 0.0848 |
| 0.75 | 174 | 0.0769 |
| 0.80 | 149 | 0.0513 |
| 0.85 | 138 | 0.0414 |
| 0.90 | 114 | 0.0236 |
| 0.95 | 93 | 0.0119 |

**Gruppo `Not_Treated` (n = 48):**

| soglia | p | zero_fraction |
| --- | --- | --- |
| 0.05 | 680 | 0.5812 |
| 0.10 | 560 | 0.5065 |
| 0.20 | 414 | 0.3814 |
| 0.30 | 334 | 0.2924 |
| 0.40 | 285 | 0.2309 |
| 0.50 | 255 | 0.1932 |
| 0.60 | 205 | 0.1282 |
| 0.70 | 175 | 0.0898 |
| 0.75 | 162 | 0.0742 |
| 0.80 | 144 | 0.0545 |
| 0.85 | 131 | 0.0425 |
| 0.90 | 109 | 0.0268 |
| 0.95 | 84 | 0.0129 |

### A.2 Soglie scelte (regola della sezione A applicata)

| Caso | Gruppo | Soglia | n | p | zero_fraction |
| --- | --- | --- | --- | --- | --- |
| R1 | all | **0.80** | 96 | 150 | 0.0590 |
| R2 | all | **0.90** | 96 | 114 | 0.0297 |
| R3 | all | **0.70** | 96 | 174 | 0.0857 |
| R4 | Treated | **0.80** | 48 | 149 | 0.0513 |
| R5 | Not_Treated | **0.80** | 48 | 144 | 0.0545 |
| L3 | all | 0.05 | 96 | 691 | 0.5851 |

Casi limite (dimensioni calcolate):

| Caso | Definizione | n | p | zero_fraction |
| --- | --- | --- | --- | --- |
| L1 | 24 campioni casuali, prev ≥ 0.20 nel sottocampione | 24 | 406–439 (media 428) | ≈0.38–0.40 |
| L2 | 100 generi casuali (nessun filtro) | 96 | 100 | 0.66–0.79 |
| L4 | 100 generi più abbondanti (conteggio totale) | 96 | 100 | 0.055 |

**Deviazione dalle soglie nominali di B (dichiarata).** Le soglie nominali
suggerite in B per `all` (≥50%, ≥20%, ≥10%) darebbero `p = 245 / 415 / 555`,
tutte **fuori** dal bersaglio `p ≈ 50–150` (R1–R2) e `150–250` (R3). Queste
soglie sono state quindi sostituite con quelle di A.2 secondo la regola
pre-dichiarata. I valori nominali non sono raggiungibili perché anche la soglia
di prevalenza più alta richiesta (50%) lascia `p = 245`, cioè **sopra** il
bersaglio di R1–R2.

### A.3 Stime dei tempi in R

Protocollo esatto (default): `spiec.easi(method="mb")`, `nlambda=20`,
`lambda.min.ratio=1e-3`, `pulsar.params=list(rep.num=20, thresh=0.05, seed=17)`,
subsample `floor(0.8·n)`. Micro-calibrazione misurata con
`spiec-easi-python/tests/real_data/calibrate_r_time.R` (gruppo `all`, n=96,
subsample 76):

| p | soglia | tempo R misurato | opt_index R | tempo Python misurato | opt_index Py |
| --- | --- | --- | --- | --- | --- |
| 82 | 0.95 | 28.5 s | 3 | — | — |
| 114 | 0.90 | 38.6 s | 3 | — | — |
| 150 | 0.80 | 51.9 s | 3 | — | — |
| 174 | 0.70 | 63.0 s | 3 | 41.4 s | 2 |
| 245 | 0.50 | 107.7 s | 3 | — | — |
| 342 | 0.30 | 184.3 s | 4 | — | — |
| 415 | 0.20 (R6) | 262.6 s | 4 | 149.3 s | 3 |
| 691 | 0.05 (L3) | 422.7 s | 4 | 439.6 s | 3 |

Un fit di potenza `T ≈ a·p^b` sui punti R misurati dà `T ≈ 0.133 · p^1.20` s
fino a `p≈245`; le misure reali di `p=342` (**184.3 s**) e `p=415`
(**262.6 s**) sono **molto sopra** l'estrapolazione, a conferma che la crescita
diventa fortemente superlineare in questo range: i valori di R6/L3 vanno
**misurati**, non estrapolati.

Nota importante: il tempo **Python è più basso di R** (solver C++ `spiecasi._mb_core`
con parallelismo su subcampioni): `p=174` → 41.4 s vs R 63.0 s; `p=415` → 149.3 s
vs R 262.6 s. L'`opt_index` di Python nella calibrazione non coincide con quello
di R perché la calibrazione usa i subcampioni casuali interni (nella batteria
gli indici di R vengono iniettati, quindi l'`opt_index` è confrontabile).

Stime derivate per seed (n=96 salvo diverso, subsample `floor(0.8·n)`):

Stime derivate (n=96 salvo diverso, subsample `floor(0.8·n)`):

| Caso | p | n | subsample | stima R | stima Py |
| --- | --- | --- | --- | --- | --- |
| R1 | 150 | 96 | 76 | ≈ 52 s | ≈ 35 s |
| R2 | 114 | 96 | 76 | ≈ 39 s | ≈ 27 s |
| R3 | 174 | 96 | 76 | ≈ 63 s | 41 s (mis.) |
| R4 | 149 | 48 | 38 | ≈ 27 s | ≈ 20 s |
| R5 | 144 | 48 | 38 | ≈ 26 s | ≈ 20 s |
| R6 | 415 | 96 | 76 | 263 s (mis.) | 149 s (mis.) |
| L1 | ≈428 | 24 | 19 | ≈ 50–55 s | ≈ 35 s |
| L2 | 100 | 96 | 76 | ≈ 34 s | ≈ 25 s |
| L3 | 691 | 96 | 76 | 423 s (mis.) | 440 s (mis.) |
| L4 | 100 | 96 | 76 | ≈ 34 s | ≈ 25 s |

Nessun caso realistico/limite pianificato (R1–R5, L1–L4) supera il tetto di
**30–60 min in R**: il caso più oneroso è `L3` (p=691, 422.7 s ≈ 7 min per
seed). I bersagli di `p`
della sezione A sono quindi guidati dalla **leggibilità statistica** e non dal
tempo; `L3` verrà comunque eseguito con un controllo di tempo (se supera
~30 min, si riduce `p` con genealogia dichiarata).

### A.4 Decisioni

- **Procedure**: si eseguono R1–R5 con 10 seed ciascuno, L1–L4 come da sezione B,
  più i sintetici S1–S4 e i patologici P1–P3.
- **Il caso `p = 1135` resta escluso dai confronti di tempo** fino all'esito del
  punto F.
- I 10 seed di R per caso saranno dichiarati nello script di generazione; gli
  indici di subcampione (0-based) e l'indice `opt.index` (1-based) saranno
  esportati e iniettati in Python.

---

## Appendice B.0 — Tempi misurati e stima totale della batteria

> Aggiunta dopo l'approvazione di A, come richiesto (punti 2–3). Le soglie di
> lettura (sezione 0, D) restano invariate.

### B.0.1 Tempi misurati (un seed, protocollo di default)

Python: `spiec-easi-python/tests/real_data/calibrate_python_time.py`
(solver C++ `spiecasi._mb_core`, `n_jobs=-1`, ~5–6 core usati).
R: `spiec-easi-python/tests/real_data/calibrate_r_time.R`.

| Caso | p | n | zero_frac | tempo R (1 seed) | tempo Python (1 seed) |
| --- | --- | --- | --- | --- | --- |
| R3 | 174 | 96 | 0.086 | 63.0 s | 41.4 s |
| **R6** | **415** | **96** | **0.380** | **262.6 s** | **149.3 s** |
| L3 | 691 | 96 | 0.585 | 422.7 s | 439.6 s |

**R6**: 262.6 s/seed ≪ tetto 30–60 min → il caso **si esegue** (nessuna riduzione
di parametri). Il caso è non degenere (`opt.index=4`, λ interno).

**L3**: 422.7 s (R) e 439.6 s (Python) per seed, entrambi ≪ tetto; anche L3 è
non degenere (`opt.index=4`).

**Python è più veloce di R** per R1–R5/R6 (≈1.6–1.8×, solver C++ con
parallelismo sui subcampioni); su L3, `p>n`, i due stack sono comparabili
(422.7 s vs 439.6 s).

### B.0.2 Stima totale della batteria (10 seed per caso)

Tempo per seed = R + Python (pipeline). Le stime non misurate usano la scalatura
osservata (`p^1.6` per R, `p^1.5` per Python; fattore `0.5` per n=48, `0.25` per
n=24). Fattore di overhead ×1.8 per i passaggi aggiuntivi del punto C (solver
isolato sul path, isolamento per subcampione, diagnostiche KKT): applicato al
totale.

| Caso | ~s/seed (R+Py) | ×10 seed |
| --- | --- | --- |
| R1 | 87 | 14.5 min |
| R2 | 66 | 11 min |
| R3 | 104 | 17.3 min |
| R4 | 47 | 7.8 min |
| R5 | 46 | 7.7 min |
| **R6** | 413 | **68.8 min** |
| L1 | 90 | 15 min |
| L2 | 59 | 9.8 min |
| **L3** | ~862 | **~144 min** |
| L4 | 59 | 9.8 min |
| S1–S4 | ~30 cad. | ~20 min (4 casi) |
| P1–P3 | ~20 cad. | ~10 min (3 casi) |

**Totale batteria (pipeline, senza overhead):** ≈ 336 min ≈ **5.6 h**.
**Con fattore di overhead ×1.8:** ≈ **10.1 h**. I due casi dominanti sono **L3**
(~2.4 h senza overhead, ~4.3 h con overhead) e **R6** (~1.15 h senza overhead,
~2.1 h con overhead).

Nessun caso supera il tetto di **30–60 min/seed in R**, quindi non scatta la
regola di stop. Il costo è però elevato: l'intera batteria a 10 seed per caso è
stimata in **~11 h** su questa macchina. Opzioni, se si vuole ridurre (da
confermare, **non** applicate in silenzio):

1. eseguire **R6 e L3 con 5 seed** (la stima scende a ~7 h);
2. eseguire la batteria **in background** (es. notturno) a 10 seed come da
   specifica;
3. eseguire prima **R1–R5, L1, L2, L4, S1–S4, P1–P3** (≈2.5 h) e **L3/R6 dopo**.

_(Fine Appendice B.0. Le sezioni B (esecuzione), C–F verranno eseguite nella fase
successiva; le soglie di lettura restano invariate.)_

---

## Appendice B.2 — Report intermedio, blocco (a)

> Blocco (a): `R1–R5, L1, L2, L4, S1–S4, P1–P3`, seed `1001…1010`, esecuzione
> **sequenziale** in background. Config di default. Python `n_jobs=-1` (thread),
> R `pulsar` **seriale** → i tempi non misurano l'efficienza del solver.
> Classificazione secondo la sezione D. Soglie di lettura **invariate**.

Dati grezzi: `tests/real_data/results/battery_results.csv`
(SHA-256 `177c6f13a30b281a7cf1ebdb809f16a2a5e04627f2fc8b6de87498e87ff5e39b`);
riepilogo: `tests/real_data/results/summary_block_a.csv`
(SHA-256 `fae0d24df0a9658a0c895230d833ae71d4cafae4142a4bdcc8242d38e3c12fd6`).

**Totale: 150 seed → 140 `identico`, 10 con solver divergente (L2; validità parziale), 0 `diverso`.**

| Caso | Tipo | n | λ uguale | rete identica | solver identico | λ interno | Classificazione | Wilson 95% (identico) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| R1 | real | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| R2 | real | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| R3 | real | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| R4 | real | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| R5 | real | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| L1 | limit | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| **L2** | limit | 10 | 10/10 | 10/10 | **0/10** | 10/10 | **solver divergente¹** | [0.000, 0.278] |
| L4 | limit | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| S1 | synth | 10 | 10/10 | 10/10 | 10/10 | **5/10** | identico | [0.722, 1.000] |
| S2 | synth | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| S3 | synth | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| S4 | synth | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| P1 | path | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| P2 | path | 10 | 10/10 | 10/10 | 10/10 | **1/10** | identico | [0.722, 1.000] |
| P3 | path | 10 | 10/10 | **4/10** | 10/10 | **1/10** | identico | [0.722, 1.000] |

¹ **L2**: solver divergente sul path; rete finale e λ identici; validità
**parziale** (vedi sezione dedicata sotto).

### Non vacuità

- **λ interno in 10/10 seed**: R1–R5, L1, L2, L4, S2, S3, S4, P1.
- **S1**: λ interno solo in **5/10** seed (i restanti 5 scelgono `λ.max`).
- **P2**: λ interno solo in **1/10** seed (9/10 scelgono `λ.max`).
- **P3**: λ interno solo in **1/10** seed (9/10 scelgono `λ.max`).

Dichiarazione esplicita: per **S1, P2, P3** il λ interno è presente solo in
**5/10, 1/10, 1/10** seed rispettivamente; nei restanti seed la scelta cade su
`λ.max` (modello vuoto/quasi vuoto) e quindi **il confronto sulla rete finale è
poco informativo** (una rete vuota coincide banalmente). Questi tre casi **non**
supportano conclusioni sulla selezione della rete.

Per **P3**, l'esito `identico` **dipende dall'eccezione `λ.max`**: in **6/10**
seed la rete finale differisce di **esattamente 1 lato** (la coppia `argmax|cor|`)
a `λ.max`, ammessa dall'eccezione documentata (`AUDIT.md` §2.E); senza
l'eccezione quei 6 seed sarebbero divergenti.

### L2 — divergenza del solver (primo passo divergente)

L2 (100 generi casuali, molto sparsi) è l'unico caso non "identico": il **solver
isola-to sul path diverge** mentre **rete finale e λ scelto coincidono** in 10/10
seed (`equivalente`). Primo passo divergente e dettaglio:

| seed | primo passo | lati diversi max sul path | `stars` raw/cum max\|diff\| |
| --- | --- | --- | --- |
| 1001 | 4_solver(sub=0, λ idx 1) | 564 | 1.23e-03 |
| 1002 | 4_solver(sub=0, λ idx 0) | 727 | 3.67e-03 |
| 1003 | 4_solver(sub=0, λ idx 2) | 883 | 2.38e-03 |
| 1004 | 4_solver(sub=0, λ idx 1) | 1160 | 1.97e-03 |
| 1005 | 4_solver(sub=0, λ idx 2) | 617 | 1.23e-03 |
| 1006 | 4_solver(sub=0, λ idx 0) | 761 | 1.38e-03 |
| 1007 | 4_solver(sub=0, λ idx 0) | 948 | 1.54e-03 |
| 1008 | 4_solver(sub=1, λ idx 2) | 895 | 1.21e-03 |
| 1009 | 4_solver(sub=0, λ idx 0) | 584 | 6.47e-04 |
| 1010 | 4_solver(sub=0, λ idx 0) | 747 | 3.58e-03 |

(λ idx = indice 0-based nella griglia; `0` = `λ.max`.) La divergenza inizia ai λ
**più regolarizzati** e prosegue a tutti i λ.

### L2b — certificazione di validità (KKT / obiettivo) e diagnostica "S di R"

Applicando il criterio dichiarato (**KKT ≤ `1e-3`**, **scarto obiettivo**
**≤ `1e-5`**) a **entrambe** le soluzioni (R e Python) usando la **S di R**, ed
eseguendo la diagnostica **"S di R"** (solver C++ alimentato con la S di R):

| Metrica | Risultato |
| --- | --- |
| divergenze totali (subcampione, λ) | **3596** |
| a `λ.max` (idx 0) | 73 — di cui **23** con ≤1 lato (eccezione), **50** con >1 lato |
| a λ **interni** | **3523** |
| **al di là dell'eccezione** (interni + `λ.max` con >1 lato) | **3573** |
| KKT ≤ `1e-3` per **entrambe** | **3596 / 3596** (max `3.89e-4`) |
| scarto obiettivo ≤ `1e-5` | **3485 / 3596** |
| **entrambe valide** (KKT ∧ obiettivo) | **3485 / 3596** |
| **NON** passano il criterio | **111 / 3596** (tutte a λ interni profondi, 1-based 14–20; **solo** per scarto obiettivo `1.00e-5…3.90e-5`; KKT regolare) |
| "S di R" fa scomparire la divergenza | **201 / 3596** (⇒ **S non è la causa**) |

**Esito L2 (certificato):** λ scelto e rete finale **identici** in 10/10 seed; il
percorso del solver è **divergente**, con soluzioni **entrambe KKT-valide**, ma
una minoranza (**111 / 3596 ≈ 3.1 %**) **non rispetta il criterio di scarto
obiettivo** dichiarato (λ profondi, scarto `≤ 3.9e-5`). La diagnostica "S di R"
esclude S come causa (recupera solo 201/3596), coerentemente con la
**degenerazione/piatteria dell'ottimo Lasso** a λ profondo (reti quasi complete).
L2 **non** è quindi "equivalente per validità" in senso pieno: va letto come
**"solver divergente, validità parziale"** (rete e λ comunque identici).

**Segno dello scarto obiettivo.** Sulle **3596** divergenze `obj_R − obj_Py` è
**simmetrico**: R migliore **1072**, Python migliore **1078**, uguali **1446**
(media `+1.2e-7`). Sulle **111** che superano la soglia: R migliore **44**,
Python migliore **67**, uguali **0** (media `+3.2e-6`, range `−3.4e-5…+3.9e-5`).
Non c'è quindi un **vantaggio sistematico** di una delle due implementazioni ⇒ è
coerente con la **piattezza/degenerazione dell'ottimo** (soluzioni multiple
quasi-ottime), **non** con un criterio di arresto o un warm-start che favorisca
uno dei due stack. La soglia `1e-5` **non** è stata modificata.

### Caveats

- `n=10` seed per caso: gli intervalli di Wilson sono ampi e **non sono una
garanzia** di equivalenza generale; "10/10" significa solo 10 successi su 10 in
quel campione (limite inferiore 95% ≈ 0.72).
- Nessun caso reale (R1–R5) è `diverso` né in errore: la regola di stop non è
  scattata.
- Le conclusioni valgono **solo** per questi casi e questa macchina.
- "Python più veloce" **non** è un confronto di efficienza (thread vs seriale).

_(Fine Appendice B.2.)_

---

## Appendice B.3 — Blocchi (b) R6 e (c) L3 (5 seed): esito ed estensione

Eseguiti in sequenza (5 seed ciascuno, `1001…1005`):

| Caso | n | λ interno | λ uguale | rete identica | solver identico | `stars` max\|diff\| | coef max\|diff\| | Classe | Wilson 95% |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| R6 (p=415) | 5 | 5/5 (idx 4) | 5/5 | 5/5 | 5/5 | 1.2e-13 | 3.7e-13 | `identico` | [0.566, 1.000] |
| L3 (p=691) | 5 | 5/5 (idx 4) | 5/5 | 5/5 | 5/5 | 2.6e-13 | 3.0e-13 | `identico` | [0.566, 1.000] |

Rete finale costante tra i seed (R6: 1617 lati; L3: 2978 lati). **Nessuna
divergenza del solver** in entrambi i casi → il criterio KKT e la diagnostica
"S di R" **non hanno oggetto** qui (non c'è nulla da certificare).

**Ipotesi (da verificare, non conclusione):** la divergenza di L2 è specifica
della presenza di generi a **prevalenza estremamente bassa (quasi-costanti)**, e
**non** della sparsità in sé. Elementi a supporto (misurati dai `counts` dei
golden):

| Caso | p | prev < 5% | prev < 10% | prev < 20% | zero_frac | min prev |
| --- | --- | --- | --- | --- | --- | --- |
| L2 | 100 | **31–47** (media ≈ 39) | 44–60 | 58–70 | 0.71–0.78 | 0.0104 |
| L3 | 691 | **0** | 136 | 276 | 0.585 | 0.052 |
| R6 | 415 | 0 | 0 | 0 | 0.380 | 0.208 |

L3 è **più sparso di R6** (`zero_frac` 0.585) e ha **136 generi con prevalenza
5–10%**, eppure è **pienamente `identico`**; L2 ha **~39 generi con prevalenza
< 5%** (min `0.0104` ≈ presente in **un solo campione**) e diverge. Coerente con
l'ipotesi, ma con `n` piccolo e **soli 2 casi** di confronto **non** è una prova:
resta **ipotesi**.

### Regola di estensione pre-registrata — applicata

L'esito dei 5 seed è **uniforme** (stessa classificazione, stesso `opt.index`,
nessun errore/divergenza) ⇒ **non ambiguo** ⇒ **si resta a 5 seed**, come da
regola a priori. Si dichiara che il campione è **piccolo**: l'intervallo di
Wilson 95% su `n=5` per "5/5 identico" è **[0.566, 1.000]** e **non è una
garanzia** di equivalenza generale.

_(Fine Appendice B.3.)_

---

## Appendice E — Riepilogo finale della batteria

Dati grezzi e riepiloghi (con SHA-256):

- `tests/real_data/results/battery_results.csv`
  (SHA-256 `d6211f91f14e6dfe402916ab47983354751614e232d6e733658b25019758ba7d`) — 160 righe.
- `tests/real_data/results/battery_summary.csv`
  (SHA-256 `df9b5aded4e5d1ad11f8a41fcc81a2407651466516666a9d511d61ebeb09f649`).
- `tests/real_data/results/validity_L2.csv`
  (SHA-256 `1acb2a73bf170ec23439ddc0ddc29d0d3e58bda05fa23b75831a4f4a152ca8eb`) — 3596 righe.

### E.1 Tabella unica per caso

| Caso | Tipo | n | λ uguale | rete identica | solver identico | λ interno | Classe | Wilson 95% (identico) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| R1 | real | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| R2 | real | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| R3 | real | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| R4 | real | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| R5 | real | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| R6 | real | 5 | 5/5 | 5/5 | 5/5 | 5/5 | identico | [0.566, 1.000] |
| L1 | limit | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| **L2** | limit | 10 | 10/10 | 10/10 | **0/10** | 10/10 | **solver divergente / validità parziale** | [0.000, 0.278] |
| L3 | real | 5 | 5/5 | 5/5 | 5/5 | 5/5 | identico | [0.566, 1.000] |
| L4 | limit | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| S1 | synth | 10 | 10/10 | 10/10 | 10/10 | **5/10** | identico | [0.722, 1.000] |
| S2 | synth | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| S3 | synth | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| S4 | synth | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| P1 | path | 10 | 10/10 | 10/10 | 10/10 | 10/10 | identico | [0.722, 1.000] |
| P2 | path | 10 | 10/10 | 10/10 | 10/10 | **1/10** | identico | [0.722, 1.000] |
| P3 | path | 10 | 10/10 | **4/10** | 10/10 | **1/10** | identico¹ | [0.722, 1.000] |

¹ P3: `identico` dipende dall'eccezione `λ.max` (6/10 seed con 1 lato = `argmax|cor|`).

**Totale: 160 seed → 150 `identico`, 10 con solver divergente (L2), 0 `diverso`.**

### E.2 Limiti dei risultati

- **Campione piccolo**: `n=10` seed (5 per R6/L3); Wilson ampi e **non** garanzia.
- **Una sola macchina / piattaforma**: Linux x86-64 (Ryzen 5 3550H); la parità è
  dimostrata solo qui.
- **Config fissa**: default (`nlambda=20`, ratio `1e-3`, `rep.num=20`, soglia
  `0.05`, regola `or`, `max_iter=10000`); non generalizzabile ad altre config.
- **Subcampioni iniettati**: gli indici di R sono iniettati in Python ⇒ si
  verifica la conseguenza deterministica a parità di subcampioni, **non** la
  riproducibilità dell'RNG di R.
- **Tempi**: Python `n_jobs=-1` (thread) vs R `pulsar` seriale ⇒ **non** misurano
  l'efficienza del solver.
- **L2**: rete e λ coincidono, ma il percorso del solver diverge e **111/3596**
  divergenze (λ profondi) superano il criterio di scarto obiettivo ⇒ **validità
  parziale**, non equivalenza piena.
- **S1/P2/P3**: parzialmente degeneri (λ interno 5/10, 1/10, 1/10) ⇒ confronto
  sulla rete poco informativo.

### E.3 Dati non verificabili

- **Altri sistemi operativi / CPU / BLAS**: non testati (scope Linux x86-64).
- **"Vecchio caso D"** citato per L2: non identificato con certezza tra i 10 seed.
- **Confronto di tempo a `p=1135`**: non eseguibile finché il lavoro non è
  dimostrato identico (punto F); **non** usato.
- **RNG di R** (StARS) non riprodotto: coperto solo via iniezione degli indici.
- **Validità KKT/obiettivo** applicata solo dove il solver diverge (L2); gli altri
  casi sono già entro `1e-10` sui passi deterministici.

### E.4 Incongruenze registrate

1. Le **soglie nominali di B** (≥50%/≥20%/≥10%) danno `p = 245/415/555`, **fuori**
   dal bersaglio ⇒ sostituite (Appendice A.2).
2. **L2**: la classificazione per-seed del CSV è `equivalente`, ma dopo la
   certificazione è **"solver divergente, validità parziale"** (111/3596 oltre il
   criterio obiettivo). Il CSV grezzo è lasciato com'è; l'etichetta corretta è nel
   testo.
3. Il **per-seed** di L2 riporta un solo "primo passo divergente" (sub 0, λ
   superficiale); la scansione completa trova **~360 divergenze/seed** su tutti i
   subcampioni e tutti i λ. Non è una contraddizione (il per-seed registra la
   prima occorrenza), ma va letto insieme alla Appendice B.2b.

_(Fine Appendice E.)_
