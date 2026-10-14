# Riepilogo — trasparenza: dati non verificabili e incongruenze

Documento di supporto a `docs/RIEPILOGO_SEMPLICE.md`. Serve a dire che cosa **non** è dimostrato
da ciò che c'è nel repository e quali contraddizioni erano presenti nei documenti, con lo stato di
ciascuna dopo il controllo del **2026-10-14**.

---

## 1. Dati che non sono verificabili (o lo sono solo in parte)

| Dato | Perché non è verificabile | Che cosa lo dimostrerebbe |
|---|---|---|
| Esito dei job CI su Windows e macOS | Nel repository c'è il file `.github/workflows/ci.yml`, che dichiara *cosa* eseguiranno i job, ma nessun esito registrato. Perciò l'equivalenza su quelle piattaforme è una promessa, non una prova | I log dei job, oppure un artefatto di esecuzione committato |
| Misure "FASE 1" (per esempio 2425.5 s, 4933.3 s, 8726.93 s, 10.740 s) | Sono registrate nei documenti e datate, ma non rieseguite ora | Riesecuzione, o gli artefatti JSON in `benchmarks/` (presenti per alcune, non per tutte) |
| Risultato delle suite lente di oggi (51 test in 1513.8 s) | È un'esecuzione locale, non un artefatto salvato; i documenti citano esecuzioni precedenti leggermente diverse (`50 passed in 1445.48 s`, `12 passed in 246 s`) | Un report committato o un job CI |
| Contenuto dei file di riferimento **prima** del 2026-10-14 | Il manifest SHA-256 dimostra la non-modifica **solo da quando esiste**. Per la storia precedente serve `git` | `git diff` contro un tag (vedi §3) |
| Idoneità degli installatori pubblicati (Windows, Linux, AppImage) | Il README li descrive; nel repository non c'è un test che li validi | Un test end-to-end sull'artefatto di release |
| Rapporti con R citati senza build di riferimento (per esempio "1.57×") | Numeri corretti ma riferiti a build diverse dello stesso caso | Citare sempre build + data (fatto il 2026-10-14 in tutti i documenti) |

Nota: dal **2026-10-14** i checksum coprono tutti e quattro gli alberi golden (`golden`,
`golden_app`, `golden_stats`, `golden_extended`), quindi la non-modifica dei file di riferimento è
verificata automaticamente **da quella data in poi**, non prima.

---

## 2. Incongruenze trovate il 2026-10-14 e stato

| # | Incongruenza | Dove | Stato |
|---|---|---|---|
| 1 | Layout di S "row-major": dato per **adottato** in §7.1 punto 6, **non adottato** in §7.6 | `docs/AUDIT.md` | **Corretta**: §7.6 ora dice che è stato adottato (commit `e73ec82`, 2026-10-13, gate su 248 problemi bit-identici) e che solo `-O3` è stato scartato |
| 2 | Costo dell'eseguibile impacchettato: "+5.3 s" rispetto al sorgente | `docs/PACKAGING.md` §5 | **Allineata**: ora riporta la misura corrente (1.504 s frozen vs 1.734 s da sorgente = −0.231 s) e dichiara il +5.3 s come numero FASE 1 superato |
| 3 | Backend dell'eseguibile impacchettato: "threads con `n_jobs=1`" | `docs/PACKAGING.md` §3, `CHANGELOG.md` (§ Fixed) | **Allineata al codice**: `python/analysis.py` usa `prefer='threads'` e `n_jobs = spieceasiJobs` se presente, altrimenti `default_n_jobs()` (core **fisici**); non esiste alcun ramo specifico per l'eseguibile impacchettato. Il test frozen ora legge il default dal log e lo confronta con `default_n_jobs()`, oltre a verificare l'uguaglianza del risultato con una corsa forzata a `n_jobs=1` |
| 4 | Default "vecchi" (nλ 25, ratio 1e-2, 50 subcampioni, frazione 0.8) elencati come divergenze, mentre §1 e §2.C li danno per allineati a R | `docs/AUDIT.md` §3 | **Corretta**: la tabella §3 ora ha una colonna "Stato" che dice per ogni voce quando è stata chiusa (Fase 3 / Fase C) o se è una scelta intenzionale ancora attiva |
| 5 | Data di aggiornamento ferma al 2026-10-13 pur contenendo la §7.7 (scritta il 14) | `docs/AUDIT.md`, intestazione | **Corretta**: intestazione aggiornata al 2026-10-14, con rimando a §7.7 |
| 6 | Rapporto con R a p=100 citato come 1.57× e come 1.54× | `docs/AUDIT.md` §7.6/§7.4, `docs/PERFORMANCE.md`, `README.md` | **Allineata**: ogni cifra ora ha data, build e condizioni (1.57× = prima del layout, 22.056 s; 1.52× = seriale dopo il layout, 21.348 s; 1.54× = `n_jobs=1` della misura con thread, 21.619 s) |

Nessuno di questi interventi ha cambiato un risultato numerico, un golden, una tolleranza o
l'eccezione stretta a λ = λ.max: hanno corretto solo lo stato dichiarato dai documenti.

---

## 3. Perché i manifest dei golden sono affidabili (e da quando)

I tre alberi `tests/golden/`, `tests/golden_app/` e `tests/golden_stats/` sono stati confrontati con
il tag **`post-cpp-extension`** (`92a18b6`) il 2026-10-14:

```
git diff --stat post-cpp-extension HEAD -- tests/golden tests/golden_app tests/golden_stats \
    ':(exclude)*/CHECKSUMS.sha256'
(nessuna differenza)
```

Quindi: **i file di riferimento di quei tre alberi sono identici** a quelli del tag; l'unica
differenza introdotta il 2026-10-14 è il file `CHECKSUMS.sha256` aggiunto in ciascun albero (che il
comando sopra esclude).

Per `golden_extended/` il confronto **non** è vuoto: 93 file, cioè le 46 matrici di path dei casi,
convertite da `.tsv` a `.tsv.gz` (46 aggiunte + 46 rimozioni; il contenuto decompresso è identico ed
è quello che registrano i checksum) più il manifest rigenerato. Per quell'albero la non-modifica del
contenuto logico è dimostrata dai checksum creati il 2026-10-13, non da `git`.

Attenzione a non invertire il senso della prova:

* il manifest certifica che il contenuto **non cambia da adesso in poi**;
* la non-modifica **rispetto al passato** è dimostrata da `git diff` contro il tag.

Le due cose insieme coprono lo stesso intervallo, ma sono due prove diverse.
