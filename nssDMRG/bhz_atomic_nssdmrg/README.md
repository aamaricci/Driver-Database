# BHZ atomico tramite nssDMRG

Il driver usa `electron_site(Hloc)` per costruire la Hamiltoniana locale
Hubbard-Kanamori. Estrae `H` dal dizionario del sito e diagonalizza con
SciFortran il settore canonico richiesto (default N=2, dimensione 6).
Non esegue sweep DMRG e non richiede checkpoint.

## Compilazione ed esecuzione

Usare `pkg-config dmrg` con la propria installazione nssDMRG e SciFortran:

```bash
export PKG_CONFIG_PATH=/Users/amaricci/opt/dmrg/gnu/cmplx/6.3.0/etc:/Users/amaricci/opt/scifor/gnu/4.26.0/etc:$PKG_CONFIG_PATH
make
./atomic_bhz FINPUT=inputATOM.conf > atomic_run.log
```

Il default e' `FC=mpif90 ARITH=cmplx`. Per la libreria reale selezionare
il suo `dmrg.pc` e ricompilare con `make -B ARITH=dble`.
`ARITH` deve corrispondere alla libreria: i tipi delle matrici locali
sono definiti dalla compilazione di nssDMRG. L'eseguibile si lancia come
processo singolo, anche quando collegato a una libreria compilata MPI.

Si puo' lavorare in una cartella separata, ad esempio:

```bash
mkdir -p run_dd
cd run_dd
cp ../inputATOM.conf .
../atomic_bhz FINPUT=inputATOM.conf JRATIO=.25 JX_FACTOR=0 JP_FACTOR=0 > atomic_run.log
```

`electron_site` stampa Hloc su standard output: per scansioni grandi
conviene redirigerlo come sopra. Il log interno del sito usa invece un
file temporaneo riscritto a ogni punto e rimosso alla chiusura.

## Input e convenzioni

Si mantengono i parametri e gli output del solver precedente.

- `NELEC`: numero di particelle del settore canonico, da 0 a 4.
- `UMIN, UMAX, NU`, `MMIN, MMAX, NM`: griglia uniforme, estremi inclusi.
  Con un solo punto si usa l'estremo inferiore.
- `JH`: Hund assoluto; `JRATIO<0` mantiene JH fissato.
  Con `JRATIO>=0`, J=JRATIO*U sostituisce JH.
- `JX_FACTOR, JP_FACTOR`: Jx/J e Jp/J. Zero per density-density;
  entrambi uno per Kanamori completa.
- `HFMODE`: passato direttamente a nssDMRG; cambia lo zero delle energie
  a numero di particelle fissato, senza alterare autostati o gap canonici.

Si usa Hloc=diag(M,-M,M,-M), Uloc=[U,U], Ust=U-2J. I coefficienti
Ust, JH, JX, JP e lo shift HFMODE sono gestiti dalla libreria.
Non si impone U>=3J; U<3J implica un coefficiente interorbitale parallelo
negativo. Non si minimizza tra diversi settori di carica.

## Operatori e proiezione

Il driver estrae i quattro operatori C dal `dot`, usando le chiavi
generate da `dot%name()` e `dot%okey`. Costruisce n=Cdag*C, S^2, Sz,
Tz=(n1-n2)/2 e D=sum_a n_a,up*n_a,down nella base completa.
Seleziona poi gli indici di N=NELEC dalla diagonale di sum(n).
Nessuna ipotesi sull'ordinamento binario della base e nessuna
ricostruzione manuale dei segni fermionici sono necessarie.
Gli operatori C non vanno proiettati prima di comporre gli osservabili,
perche' cambiano il numero di particelle.

Il problema atomico attuale ha matrici reali anche con la libreria
complessa; il driver controlla che le parti immaginarie siano nulle
prima di convertirle e usare `eigh` reale. Controlla anche hermiticita',
conservazione di N e diagonalita' dell'operatore numero.

L'inizializzazione configura solo lo stato di libreria richiesto da
`electron_site` (due orbitali, due spin, modo normal, PBC disattivate,
interazioni, HFMODE e log). Non usa `read_input` del DMRG, per evitare
di introdurre parametri di catena e sweep in un calcolo di singolo sito.

## Output e verifiche

`atomic_ground.dat` contiene U, M, J, E0, degenerazione, gap sopra il
sottospazio fondamentale, S^2, Sz, Tz, D, n1, n2.
Gli osservabili sono medie equiprobabili sul sottospazio fondamentale,
con tolleranza 1e-10*max(1,max|E|). Se l'intero settore e' degenere il
gap e' scritto come zero. Il gap e' canonico, non un gap di carica.
`atomic_spectrum.dat` contiene U, M, J, indice del livello, energia.
Entrambi vengono sovrascritti nella directory di esecuzione.

Prima dello shift HFMODE, i livelli N=2 sono:

```
U-3J (due volte), U-2J-Jx, U-2J+Jx,
U-sqrt(4M^2+Jp^2), U+sqrt(4M^2+Jp^2).
```

All'avvio il driver verifica questi livelli e lo spin del fondamentale
Hund in due punti di riferimento, con e senza SF/PH.
Per controllare spettri analitici, degenerazioni e osservabili su tutti
i settori, entrambi gli HFMODE e i due modi di fissare Hund:

```bash
python3 check_atomic.py
# Confronto facoltativo con un precedente eseguibile:
python3 check_atomic.py --reference /percorso/al/solver/precedente
```

Il controllo Python richiede NumPy. Le formule analitiche sono usate
solo nei controlli; la Hamiltoniana del calcolo viene dal `dot`.

Per J/U=0.25 senza SF/PH, la frontiera a N=2 e' |M|=3U/8.
I dati e il PNG gia' presenti in `phase-diagram/` restano utilizzabili
con gli stessi formati di output.
