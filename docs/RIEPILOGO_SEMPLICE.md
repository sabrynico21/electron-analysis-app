# Riepilogo semplice — che cosa è stato verificato, perché, e cosa manca

Questo documento riassume, in linguaggio semplice, i controlli fatti sulla parte "statistica"
del progetto. Ogni numero viene da un file del repository (`docs/AUDIT.md`,
`docs/PERFORMANCE.md`, `docs/STATISTICAL_EQUIVALENCE.md`, `CHANGELOG.md`, `docs/PACKAGING.md`,
`README.md`) o da un test eseguito su questa macchina. I dati non verificabili e le
incongruenze trovate nei documenti sono in **`docs/RIEPILOGO_TRASPARENZA.md`**.

Il progetto è un'applicazione desktop che, da una tabella di abbondanza di microrganismi,
calcola una rete di co-occorrenza. Uno dei metodi disponibili è **SPIEC-EASI (metodo MB)**.
Il repository contiene anche una libreria Python separata (`spiec-easi-python`) che è una
traduzione di quel metodo, scritto in origine nel linguaggio R.

---

## 1. Che cosa significa "fedele alla versione R"

- **R** è il programma di riferimento (un linguaggio usato in statistica). Il riferimento è
  composto da tre pacchetti: **SpiecEasi 1.99.0**, **huge 2.0.1** e **pulsar 0.3.13**.
- "Fedele" significa: **a parità di dati in ingresso, la versione Python esegue gli stessi
  passaggi, nello stesso ordine, e arriva allo stesso risultato**. Il codice Python è riscritto
  da zero, quindi non è una copia del codice R.
- Come si controlla: confrontando **passo per passo**, non solo il risultato finale. R salva i
  risultati intermedi in file chiamati **golden** ("di riferimento") con 17 cifre; i test
  verificano che Python riproduca quegli stessi numeri. Se un passo non torna, il test si ferma
  e dice **quale** passo ha fallito e di quanto.
- Un punto delicato: la scelta casuale dei "subcampioni" (StARS, punto 5 della tabella) non è
  riproducibile in Python, perché il generatore di numeri casuali di R è diverso. Nei confronti
  i subcampioni di R vengono quindi **impostati a mano** dentro Python: così si isolano tutti i
  passaggi deterministici. L'unico test che **non** usa questo trucco è quello di equivalenza
  statistica (punto 13 della tabella).
- Le **tolleranze** ("quanto scarto accettiamo") sono dichiarate **prima** di eseguire i test
  e non sono mai state allargate per far passare un test. Nei file golden non è stato cambiato
  nulla.

> Nota sulla numerazione: la richiesta salta dal punto 1 al punto 3 (il punto 2 non c'è).
> Ho tenuto la numerazione della richiesta, quindi la tabella è il punto 3 e i punti
> successivi seguono il testo originale.

---

## 3. Tabella dei test

Legenda esito: "passa" = esito dell'ultima esecuzione su questa macchina
(153 test in totale: 102 veloci + 51 lenti).

| Cosa è stato testato | Perché è stato fatto | Come (in una riga) | Risultato breve |
|---|---|---|---|
| 1. Trasformazione dei dati (CLR) | Primo passo: un errore si trascina su tutto | Confronto dei valori trasformati sugli stessi dati | Passa. Scarto ≤ 2e-15 (fase 2), 1.3e-15…3.6e-15 sui 23 casi estesi |
| 2. Parametri da provare (λ massimo e griglia) | Da qui dipende la ricerca della rete | Confronto del valore massimo e della griglia generati | Passa. Scarto 0…3.3e-16 (λ massimo) e 0…1.7e-16 (griglia) |
| 3. Solver a valore fisso (i coefficienti) | È il cuore del calcolo | Confronto dei coefficienti a un valore fissato | Passa. Scarto ≤ 1.1e-13 (fase 2), 6.7e-15…7.4e-13 sui casi estesi |
| 4. Simmetrizzazione e pesi dei lati | La rete finale non è orientata | Confronto dei lati e dei pesi con la regola di R | Passa. 0 lati di differenza a valore fisso; pesi = media dei due versi |
| 5. Scelta del parametro finale (StARS) | Decide "quanto fitta" la rete | Confronto della curva di instabilità e dell'indice scelto | Passa. Curva ≤ 1e-15 (grezza e cumulata); indice e λ sempre uguali |
| 6. Casi di confine: soglia eguagliata, λ massimo | Gli errori `<` invece di `≤` si vedono solo ai bordi | 6 casi costruiti con R + il caso λ massimo | Passa (6/6). A λ massimo al più 1 lato di differenza, e solo la coppia più correlata |
| 7. Dataset estesi: 23 casi (seed, 75% zeri, più specie che campioni, n > 144) | Verifica fuori dai casi "comodi" | Confronto degli 9 passi, caso per caso | Passa. Passi 1–8 identici; al passo 9 **1 lato diverso in 8 casi su 20**, tutti a λ massimo |
| 8. Casi grandi (p = 100) | Dimensione realistica | 3 casi eseguiti per intero | Passa (3 casi). **0 differenze** in tutti i passi |
| 9. App: colonne costanti, regola per n > 144 | L'app deve comportarsi come R sui dati reali | Funzione dell'app confrontata con R + controllo dei parametri passati | Passa (8 test). Colonna costante scartata (R darebbe errore); regola per n > 144 verificata |
| 10. Casi mal condizionati + criterio alternativo | Lì la soluzione **non è unica** | Verifica di "validità" con due soglie dichiarate prima | Passa. 2.478e-04 (limite 1e-3) e 2.390e-06 (limite 1e-5) |
| 11. Velocità e determinismo | Con più thread il risultato non deve cambiare | Stesso caso con 1/2/4/8 lavoratori, confronto esatto | Passa. **Identici bit per bit** (anche col metodo a processi) |
| 12. Eseguibile impacchettato | È ciò che scarica l'utente | Si lancia l'eseguibile e si leggono backend e lavoratori dal log | Passa (2 test). Nucleo compilato, **thread + core fisici**, risultato = corsa a 1 lavoratore |
| 13. Equivalenza statistica senza gli indici di R | Ognuno estrae i propri subcampioni (uso reale) | 2 casi × 30 seed: curve medie e reti ottenute | Passa. Curve entro 0.011 e 0.0096; reti vuote (proprietà dei dati) |
| 14. Checksum dei file di riferimento | Nessun cambiamento nascosto | SHA-256 ricalcolato e confrontato con la lista salvata | Passa. 464 voci su 4 alberi, tutte combaciano (**da oggi in poi**) |

Totale: **153 test** (102 veloci + 51 lenti). Ultima esecuzione completa sullo stato attuale:
102 veloci in 21.0 s e 51 lenti in **1414.1 s (23 min 34 s)**, tutti passati.

---

## 4. Cosa abbiamo scoperto e corretto lungo il percorso

1. **Il solver veniva usato "a freddo" invece che "a caldo".** R risolve tutti i valori della
   griglia in una sola passata e riusa il risultato precedente come punto di partenza; Python
   ricominciava da zero ogni volta. In casi difficili si arrivava a reti diverse (fino a 4e-2).
   Corretto: ora Python fa come R e quei casi combaciano.
2. **Alcuni valori di default non erano quelli di R.** La libreria usava 25 valori di λ invece
   di 20, un intervallo 1e-2 invece di 1e-3 e 50 subcampioni invece di 20; l'app forzava la
   frazione di subcampione a 0.8 anche per n > 144, dove R usa un'altra regola. Allineato a R.
3. **Un confronto al bordo era sbagliato di uno.** La regola che sceglie il valore finale
   usava `≤` invece di `<`: quando la curva toccava **esattamente** la soglia, il risultato
   cambiava. Corretto, con 6 casi generati con R.
4. **L'eseguibile impacchettato si rompeva a ogni analisi SPIEC-EASI.** Il sistema dei
   "processi" riavviava l'eseguibile con un parametro che l'eseguibile rifiutava. Corretto
   usando i "thread" (che non riavviano nulla); un test di regressione lo verifica.
5. **Un test che non verificava nulla, e un confronto di tempi non equo.** Lo scenario di
   "variabilità" sceglieva sempre λ massimo, cioè il blocco indeterminato: confrontava reti
   vuote e sembrava passare. Ora ci sono tre scene con segnale forte e un controllo che
   fallisce se i casi utili sono troppo pochi. Inoltre il pacchetto R **huge** è compilato con
   OpenMP (parallelismo interno, 1.56×–2.13× sul solo solver): a parità di risorse il nostro
   nucleo è alla pari (p=100: 1.66 s contro 1.86 s a un thread).

---

## 5. Risultati sulle prestazioni (prima e dopo)

Macchina: AMD Ryzen 5 3550H, **4 core fisici / 8 logici**, Linux 6.18.53 x86_64, glibc 2.44,
Python 3.14.7, R 4.6.1. Pipeline completa (StARS + fit finale), mediana di 3 esecuzioni, un solo
lavoratore su entrambi i lati.

| caso | solo Python (prima) | Python con nucleo compilato | R (riferimento) |
|---|---|---|---|
| p=30, n=100 | 97.819 s | 0.572 s | 2.108 s |
| p=50, n=50 | 453.408 s | 3.991 s | 4.450 s |
| p=100, n=80 | 1565.844 s | 22.056 s | 14.040 s |

Dopo un'ulteriore ottimizzazione dei dati in memoria (adozione del "row-major", risultati
identici bit per bit): 0.514 s / 3.691 s / 21.348 s.

Con più lavoratori ("n_jobs"), stessa macchina:

| | p=30 | p=50 | p=100 |
|---|---|---|---|
| 1 lavoratore | 0.516 s | 3.804 s | 21.619 s |
| 4 lavoratori | 0.239 s | 1.497 s | **8.734 s** |
| 8 lavoratori | 0.196 s | 1.250 s | **6.861 s** |
| R | 2.108 s | 4.450 s | 14.040 s |

Altri numeri misurati:

- la suite di test veloce passa da 662 s a 18.66 s;
- l'avvio dell'eseguibile impacchettato su un caso piccolo passa da 10.740 s a **1.504 s**
  (di cui circa 9.2 s erano il vecchio solver e una libreria importata inutilmente);
- a p = 100 il 99.7% del tempo è dentro il nucleo del solver; tutto il resto del Python pesa
  lo 0.1%.

Avvertenza: a p=100 e a un lavoratore il nostro nucleo resta ~1.5× più lento di R, ma il
confronto non è alla pari, perché R usa OpenMP al proprio interno.

---

## 6. Cosa **non** è garantito

Va detto chiaramente, senza addolcire:

- **Non è una garanzia generale.** L'equivalenza è dimostrata **solo sui casi testati** (le
  dimensioni, i seed e i tipi di dati elencati nei documenti).
- **Solo su Linux x86-64 con compilatore GCC.** Windows e macOS sono previsti nella
  configurazione CI, ma **non sono provati**: finché quei job non girano e passano, non
  possiamo dire che si comportino allo stesso modo. Compilatori e processori diversi possono
  cambiare le ultime cifre dei numeri.
- **Senza gli indici di R, il risultato è statisticamente simile, non identico.** Quando
  entrambi i programmi estraggono i propri subcampioni, le reti possono differire; il test
  mostra che le curve di instabilità sono vicine (scarto massimo ~0.01), non uguali.
- **A λ massimo il risultato è indeterminato anche in R.** Nello stesso punto, con dati
  diversi, R a volte tiene un lato e a volte nessuno. Python in quel punto restituisce sempre
  la rete vuota, che è la scelta matematicamente corretta. Nei test è ammessa solo questa
  eccezione, ed è stretta.
- **Le differenze numeriche esistono.** Il nucleo compilato **non** è identico bit per bit né
  al riferimento Python né a R: 0 problemi su 244 sono bit-identici, e lo scarto massimo è
  9.470e-14 (cioè ~1e-13) nei casi normali e 1.652e-13 rispetto ai file di riferimento.
- **I casi "degeneri" restano diversi.** In una famiglia costruita apposta (dove la soluzione
  non è unica) lo scarto massimo è 1.881e-01 su 25 problemi e 40 blocchi di λ. R non sceglie
  né l'una né l'altra soluzione: il problema è dei dati, non del codice. Questi casi sono
  verificati con un criterio di validità, non con un confronto di numeri.
- **Deviazioni volute dell'app:** per dataset molto piccoli e con moltissime specie
  (n ≤ 10 e p ≥ 200) l'app riduce il numero massimo di iterazioni da 10000 a 3000 (R non lo
  fa), e scarta le colonne costanti dove R darebbe un errore.
- **Due funzioni di R non sono implementate:** la sostituzione "moltiplicativa" degli zeri e
  la griglia di λ lineare.
- **I checksum dimostrano la non-modifica solo da adesso in poi.** Coprono tutti e quattro gli
  alberi golden (464 voci), ma un manifest non dice nulla sulla storia precedente di un file:
  per quella serve `git` (vedi `docs/RIEPILOGO_TRASPARENZA.md` §3).

Quindi: non va scritto "identico a R" né "verificato al 100%". La frase corretta è
"equivalente sui casi testati, con i limiti elencati".

---

## 7. Cosa resta da fare (in ordine)

1. **Far girare la CI su Windows e macOS** (Intel e ARM): build e sweep completa contro R.
2. **Verifiche manuali prima di una release:** clone pulito, installazione da zero, esecuzione
   dell'eseguibile impacchettato (provato finora solo su Linux). Checksum, percorsi personali
   e sweep sono già automatici.
3. **Decidere la licenza: è bloccante.** Il repository dichiara **MIT**, ma i tre pacchetti di
   riferimento sono **GPL** (GPL-2 `huge`, GPL-3+ SpiecEasi e pulsar): essendo una traduzione
   di quei metodi, la licenza attuale non è compatibile e non si può pubblicare così.
4. **Repository privato, poi pubblico**, una volta risolta la licenza.
5. **Integrare la libreria aggiornata nell'app** e rieseguire i test di integrazione.
6. **Provare su dati reali** (finora tutti i confronti sono su dati generati).
7. **Eventuale ottimizzazione futura:** parallelismo OpenMP nel solver (R lo usa). Rimandata,
   nessun codice.

---

## 8. Sintesi in una frase

Oggi possiamo dire che la versione Python riproduce il metodo di R in modo **fedele e
verificato sui casi testati su Linux x86-64**, con differenze numeriche di ordine 1e-13, un
limite dichiarato a λ massimo (indeterminato anche in R) e una famiglia di casi degeneri
riportata senza nasconderla; restano da fare la prova su Windows e macOS, la scelta della
licenza e i test su dati reali.

---

## 9. Trasparenza

I dati non verificabili e le incongruenze trovate nei documenti — con lo stato di ciascuna dopo il
controllo del 2026-10-14 — sono in **`docs/RIEPILOGO_TRASPARENZA.md`**.

