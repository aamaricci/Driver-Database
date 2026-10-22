# Driver BHZ di post-processing

Compilazione con nssDMRG **CMPLX**, usando gli stessi flag MPI della libreria:

```sh
make structure polarization
# Per scrivere gli eseguibili in una directory locale già esistente:
make structure polarization DIREXE=.
```

Eseguire dalla directory contenente i checkpoint del ground state e le
rotazioni salvate (`save_block=T`, `save_umat=T` durante il calcolo originale):

```sh
structure_bhz_1d FINPUT=used.DMRG.conf
polarization_bhz_1d FINPUT=used.DMRG.conf
```

In MPI usare il consueto `mpirun -np ...`. Nessun driver chiama `run_DMRG`;
`IRUN` non viene usato. `MH` e `LAMBDA` vengono letti dal file, con gli stessi
default del driver BHZ (0.5 e 0.3). Servono `Nspin=Norb=2`, `QNDIM=2`,
`DMRG_MODE=normal` e parametri coerenti con il checkpoint. In nssDMRG 6.3.0
il settore usa `QN_DENSITY=1,1` e `QN_OFFSET=0,0` a metà riempimento:
i vecchi record `QNTYPE` e `DMRG_QN` non sono più accettati.
Se i checkpoint sono archiviati in `restart_*.tgz`, estrarli e predisporre
la cartella `restart/` prima di lanciare i driver.

## Fattori di struttura

`structure_bhz_1d.f90` scrive `Sq` e `Nq` con il suffisso `label_DMRG('u')`,
come gli altri driver. Colonne: q in radianti/sito, parte reale e immaginaria.

Si misura F(q)=sum_ij exp[i q(j-i)] C_ij/M, inclusi i=j, su M siti;
S(q) usa Sz=(n_up-n_down)/2, N(q) usa n=n_up+n_down.
S(q) è la componente longitudinale, senza assumere simmetria SU(2).

Input opzionale:

```text
BULK_CUT=0
STRUCTURE_CONNECTED=F
```

`BULK_CUT` elimina quel numero di siti a ciascun bordo solo per S(q), N(q).
La griglia è q=2*pi*k/M, k=0,...,M-1. In OBC è una trasformata a onde piane.
Il default conserva le correlazioni raw. `STRUCTURE_CONNECTED=T` sottrae
localmente <A_i><A_j>, rimuovendo anche eventuali contributi di ordine CDW.
Nella catena intera a carica fissata, N(0) connected deve annullarsi;
a metà riempimento N(0) raw=4L. S(0) connected si annulla a Sz totale fissato.

## Polarizzazione

`polarization_bhz_1d.f90` scrive `z_sigma` con lo stesso suffisso.
Colonne: spin (1 up, 2 down), Re z, Im z, |z|, arg z [rad], P=arg z/(2*pi).
Si calcola z_sigma=<exp[i 2*pi/L sum_j j*n_j,sigma]> sull'intera catena,
j=1,...,L, senza sottrazione di background. P è definita modulo 1;
la fase non è significativa quando |z| tende a zero.

Il marker conserva la convenzione di origine e il fattore di parità della
lunghezza: non viene arrotondato a un indice discreto. L'interpretazione
topologica richiede la convenzione appropriata; in OBC il risultato dipende
anche dallo stato dei bordi. I due spin sono mantenuti separati perché il
marker della carica totale può perdere l'informazione dei singoli settori.
