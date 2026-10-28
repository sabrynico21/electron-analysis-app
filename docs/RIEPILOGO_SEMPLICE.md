
# Riepilogo semplice — i test eseguiti e i limiti

Questo documento riassume come è stato controllato che la versione Python
del metodo **SPIEC-EASI (MB)** si comporti come la versione di riferimento scritta in R. Ogni
numero viene dai documenti del progetto (`docs/AUDIT.md`, `docs/PERFORMANCE.md`,
`docs/STATISTICAL_EQUIVALENCE.md`, `CHANGELOG.md`, `README.md`) o da un test eseguito su questa
macchina.

## Che cosa significa "fedele alla versione R"

- **R** è il programma di riferimento. Il riferimento è composto da tre pacchetti:
  **SpiecEasi 1.99.0**, **huge 2.0.1** e **pulsar 0.3.13**. Il codice Python è una riscrittura
  di quei pacchetti, non una copia.
- "Fedele" significa: **a parità di dati in ingresso, la versione Python esegue gli stessi
  passaggi, nello stesso ordine, e arriva allo stesso risultato.**
- Come si controlla: confrontando **passo per passo**, non solo il risultato finale. R salva i
  risultati intermedi in file chiamati **golden** ("di riferimento") con 17 cifre; i test
  verificano che Python riproduca quegli stessi numeri. Se un passo non torna, il test si ferma
  e dice **quale** passo ha fallito e di quanto.
- Un punto delicato: la scelta casuale dei "subcampioni" (StARS, punto 5 della tabella) non è
  riproducibile in Python, perché il generatore di numeri casuali di R è diverso. Nei confronti
  i subcampioni di R vengono quindi **impostati a mano** dentro Python: così si isolano tutti i
  passaggi deterministici. L'unico test che **non** usa questo trucco è quello di equivalenza
  statistica (punto 13 della tabella).
- Le **tolleranze** sono dichiarate **prima** di eseguire i test
  e non sono mai state allargate per far passare un test. Nei file golden non è stato cambiato
  nulla.

---

## Tabella dei test

I test sono **153** in tutto: 102 veloci e 51 lenti.

### A. Confronto con R, passo per passo

| Cosa è stato testato | Perché | Come | Esito |
|---|---|---|---|
| **1. Trasformazione dei dati (CLR)** | È il primo passo: un errore qui si trascina su tutto il resto | Si confrontano i valori trasformati sugli stessi dati | Passa. Scarto ≤ 2e-15 (fase 2); 1.3e-15…3.6e-15 sui 23 casi estesi |
| **2. Parametri da provare** (λ massimo e griglia) | Da questi valori dipende tutta la ricerca della rete | Si confrontano il valore massimo e la griglia generati dai due programmi | Passa. Scarto 0…3.3e-16 (λ massimo) e 0…1.7e-16 (griglia) |
| **3. Solver a valore fisso** (i coefficienti) | È il cuore del calcolo: qui R usa una routine in C++ e Python una riscrittura | Si confrontano i coefficienti calcolati a un valore fissato | Passa. Scarto ≤ 1.1e-13 (fase 2); 6.7e-15…7.4e-13 sui casi estesi |
| **4. Simmetrizzazione e pesi dei lati** | Il calcolo produce una tabella *orientata*, ma la rete di co-occorrenza non ha versi: i lati vanno resi simmetrici. Se la regola di combinazione cambia, cambia la rete | Si prende la tabella orientata, si applica la regola di R (lato presente se **almeno uno dei due versi** è diverso da zero; peso del lato = **media dei due versi**) e si confronta con R sia l'elenco dei lati sia i pesi | Passa. Elenco dei lati **identico** a R (0 lati di differenza) su tutti i casi confrontati; pesi pari alla media dei due versi, come in R |
| **5. Scelta del parametro finale (StARS)** | È il passo che decide "quanto fitta" la rete | Si confrontano la curva di instabilità e l'indice scelto | Passa. Curva ≤ 1e-15 (grezza e cumulata); indice e λ sempre uguali |

### B. Casi difficili e dimensioni realistiche

| Cosa è stato testato | Perché | Come | Esito |
|---|---|---|---|
| **6. Casi di confine:** soglia raggiunta esattamente, e λ massimo | Gli errori di confronto (`<` invece di `≤`) si vedono solo ai bordi | 6 casi costruiti apposta con R, più il caso λ massimo | Passa (6/6). A λ massimo è ammesso **al più 1 lato** di differenza, e solo sulla coppia più correlata |
| **7. Dataset estesi:** 23 casi (più seed, 75% di zeri, più specie che campioni, campioni > 144) | Serve a vedere se il metodo regge anche fuori dai casi "comodi" | Per ogni caso si confrontano gli stessi 9 passi | Passa. Passi 1–8 identici; al passo 9 (rete finale) **1 lato diverso in 8 casi su 20**, tutti a λ massimo |
| **8. Casi grandi (p = 100)** | È la dimensione realistica: serve la prova, non la stima | 3 casi eseguiti per intero, tutti i passi confrontati | Passa (3 casi). **0 differenze** in tutti i passi |
| **10. Casi mal condizionati** (soluzione non unica) + criterio alternativo | Lì confrontare i numeri non ha senso: servono altre domande | Si verifica che **entrambe** le soluzioni siano valide per lo stesso problema, con soglie dichiarate prima | Passa. 2.478e-04 (limite 1e-3) e 2.390e-06 (limite 1e-5) |

### C. Il programma vero: app, imballaggio, robustezza, file

| Cosa è stato testato | Perché | Come | Esito |
|---|---|---|---|
| **9. Comportamento dell'app** (colonne costanti, regola di subcampionamento per n > 144) | L'app deve comportarsi come R sui dati reali degli utenti | Si esegue la funzione dell'app e si confronta con R, più controlli sui parametri passati | Passa (8 test). La colonna costante viene scartata (R darebbe errore); la regola per n > 144 è verificata |
| **11. Velocità e determinismo** | Il calcolo usa più "thread": il risultato non deve cambiare | Stesso caso con 1, 2, 4, 8 lavoratori, confronto esatto dei risultati | Passa. **Identici bit per bit** (anche col vecchio metodo a processi) |
| **12. Eseguibile impacchettato** (quello che si scarica) | Se è rotto, l'app non funziona per nessuno | Si lancia l'eseguibile vero e si leggono dal suo log backend e numero di lavoratori | Passa (2 test). Usa il nucleo compilato, **thread + core fisici**, risultato identico a una corsa a 1 lavoratore |
| **13. Equivalenza statistica senza gli indici di R** | È la verifica più severa: ognuno estrae i propri subcampioni, come nell'uso reale | 2 casi × 30 seed: si confrontano le curve medie e le reti ottenute | Passa. Curve entro 0.011 e 0.0096; reti vuote (proprietà dei dati, non del codice) |
| **14. Checksum dei file di riferimento** | Serve a impedire che un file di riferimento cambi di nascosto | Si ricalcola l'impronta SHA-256 di ogni file e si confronta con la lista salvata | Passa. 464 voci su 4 gruppi di file, tutte combaciano (**prova valida da oggi in poi**) |

**Totale:** 153 test (102 veloci + 51 lenti). Ultima esecuzione completa sullo stato attuale:
102 veloci in 21.0 s, 51 lenti in **1414.1 s (23 min 34 s)**, tutti passati.

---

## Prestazioni: la versione attuale contro R


Macchina di misura: AMD Ryzen 5 3550H (**4 core fisici / 8 logici**), Linux x86_64, Python 3.14.7,
R 4.6.1. Pipeline completa (StARS + fit finale), mediana di 3 esecuzioni, **un solo lavoratore su
entrambi i lati** (confronto alla pari).

| caso | versione attuale (compilata) | R | confronto |
|---|---|---|---|
| p=30, n=100 | 0.514 s | 2.108 s | **4,1× più veloce di R** |
| p=50, n=50 | 3.691 s | 4.450 s | **1,2× più veloce di R** |
| p=100, n=80 | 21.348 s | 14.040 s | **1,5× più lento di R** |


**Con più lavoratori** (nome tecnico "n_jobs"; il valore predefinito è il numero di core fisici,
quindi 4 su questa macchina):

| lavoratori | p=30 | p=50 | p=100 | confronto con R a p=100 |
|---|---|---|---|---|
| 1 | 0.516 s | 3.804 s | 21.619 s | 1,5× più lento |
| **4 (predefinito)** | 0.239 s | 1.497 s | **8.734 s** | **1,6× più veloce** |
| 8 | 0.196 s | 1.250 s | **6.861 s** | 2,0× più veloce |
| R (riferimento) | 2.108 s | 4.450 s | 14.040 s | — |

Nota: le due tabelle vengono da due serie di misure diverse (stessa macchina e stesse condizioni),
e per questo il caso p=100 a un lavoratore compare come 21.348 s in una e 21.619 s nell'altra.


**Perché a p=100, con un solo lavoratore, restiamo più lenti di R.** R (pacchetto `huge`) è
compilato con OpenMP, cioè usa più core *al proprio interno*: sul solo solver, a p=100 passa da
1.862 s a un thread a 0.873 s con gli 8 thread predefiniti (circa 2,1× più veloce). Quindi il
confronto "un lavoratore contro un lavoratore" non è alla pari. A parità di risorse il nostro
nucleo è alla pari (p=100: 1.66 s contro 1.86 s a un solo thread). Usare OpenMP anche nel nostro
codice è una proposta futura, **non implementata**.

Una piccola ottimizzazione interna già adottata (la disposizione dei dati in memoria, detta
"row-major") vale circa il 6% e produce risultati **identici bit per bit**: non cambia nulla nei
numeri confrontati con R.

---

## Cosa **non** è garantito

Questa è la parte più importante del documento. In una frase: **il comportamento è verificato sui
casi testati e su questa piattaforma; non è una garanzia universale.** Qui sotto ogni limite, con
la sua conseguenza pratica.

1. **Non è una garanzia generale.** I confronti coprono le dimensioni, i seed e i tipi di dati
   elencati nei documenti. *In pratica:* su dati diversi — per esempio con molte più specie di
   quelle provate — non sappiamo come si comporta, e non si può promettere che sia uguale a R.
2. **Solo Linux x86-64 con compilatore GCC.** Windows e macOS sono previsti nella CI, ma quei job
   non sono ancora stati eseguiti. *In pratica:* fino ad allora l'affermazione "uguale a R" vale
   **solo per Linux**; su Windows e macOS non sappiamo se le ultime cifre coincidono.
3. **Senza gli indici casuali di R il risultato è statisticamente simile, non identico.** Nel
   confronto più severo ogni programma estrae i propri subcampioni, quindi le reti finali possono
   differire. *In pratica:* su dati reali due esecuzioni possono dare reti leggermente diverse
   (curve di instabilità entro circa 0.01): è normale, la procedura è casuale per costruzione.
4. **A λ massimo il risultato è indeterminato anche in R.** In quel punto preciso, con dati
   diversi, R a volte tiene un lato e a volte nessuno. *In pratica:* lì Python restituisce sempre
   la rete vuota (scelta voluta e matematicamente corretta) e nei test è ammessa **solo** questa
   differenza, **al più un lato**, sulla coppia più correlata.
5. **Le differenze numeriche esistono, piccole ma reali.** Il nucleo compilato **non** è identico
   bit per bit al riferimento Python né a R: su 244 problemi, **0** sono bit-identici; lo scarto
   massimo sui coefficienti è 9.470e-14 nei casi normali e 1.652e-13 rispetto ai file di
   riferimento. *In pratica:* "9e-14" è un errore sulla quattordicesima cifra decimale: la rete
   scoperta è la stessa, cambiano solo gli ultimi decimali dei pesi (non i lati).
6. **I casi "degeneri" restano diversi.** In una famiglia costruita apposta (dati quasi
   collineari, dove la soluzione matematica non è unica) lo scarto massimo è 1.881e-01 su 25
   problemi e 40 blocchi di λ. *In pratica:* lì la scelta della soluzione è arbitraria anche in R,
   che non privilegia né la nostra né l'altra; per questo quei casi si verificano con un criterio
   di validità ("entrambe sono soluzioni accettabili") invece che confrontando i numeri.
7. **Deviazioni volute dell'app.** Con dataset molto piccoli e molto "larghi"
   (n ≤ 10 campioni e p ≥ 200 specie) l'app riduce il numero massimo di iterazioni da 10000 a
   3000; inoltre scarta le colonne costanti, dove R darebbe un errore. *In pratica:* in quelle
   situazioni l'app è più permissiva e più veloce di R, quindi non identica.
8. **Due funzioni di R non sono implementate:** la sostituzione "moltiplicativa" degli zeri e la
   griglia di λ lineare. *In pratica:* sono opzioni che R offre e Python no; se servono, non sono
   disponibili.

---

## Sintesi in una frase

La versione in Python riproduce il metodo di R in modo **fedele e verificato sui casi testati su
Linux x86-64**; le differenze numeriche sono di ordine 1e-13 (la rete scoperta è la stessa), con
un limite dichiarato a λ massimo (indeterminato anche in R) e una famiglia di casi degeneri
riportata senza nasconderla. Restano da fare la prova su Windows e macOS, la scelta della licenza
e i test su dati reali.

---



