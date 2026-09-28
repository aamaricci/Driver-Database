#!/usr/bin/env python3
"""Extract and plot BHZ entanglement spectra from a DMRG data directory.

Usage: python3 bhz_entanglement_analysis.py --root /path/to/M1.7

The input directory is expected to contain subdirectories named U* with files
LambdaQ_left/right_L*_*.dmrg.  Each LambdaQ file stores blocks of the form

    q Nq
    (lambda_1,0)
    ...

The script writes numerical tables and publication-quality PNG figures into
ROOT/entanglement_spectrum.  The plotting backend is Agg, so no GUI is needed.
"""
import argparse, math, re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read_lambda_file(path):
    """Read one symmetry-resolved RDM spectrum.

    Returns tuples (q, index-within-q, complex-eigenvalue).  The Fortran
    complex syntax is parsed explicitly because it is not valid Python syntax.
    """
    values=[]; q=None; expected=0
    for raw in path.read_text().splitlines():
        line=raw.strip()
        if not line: continue
        if expected == 0:
            f=line.split(); 
            q,expected=int(f[0]),int(f[1]); 
            continue
        m=re.match(r"^\(\s*([^,]+)\s*,\s*([^\)]+)\s*\)$",line)
        if not m: 
            raise ValueError(f"Invalid eigenvalue in {path}: {line}")
        z=complex(float(m.group(1).replace("D","E")),float(m.group(2).replace("D","E")))
        iq=sum(x[0] == q for x in values)+1
        values.append((q,iq,z)); 
        expected-=1
    if expected: raise ValueError(f"Truncated sector q={q} in {path}")
    return values



def extract(path, side, outdir):
    """Convert lambda eigenvalues into xi=-log(lambda) and save one table."""
    rows=[]
    for q,iq,z in read_lambda_file(path):
        if abs(z.imag)>1e-9*max(1.,abs(z.real)): raise ValueError(f"Complex lambda in {path}")
        if z.real>0: rows.append((z.real,q,iq))
    rows.sort(reverse=True); 
    out=outdir/f"spectrum_{side}_{path.parent.name}.dat"
    with out.open("w") as f:
        f.write("# rank xi=-log(lambda) lambda q index_in_sector\n")
        for rank,(lam,q,iq) in enumerate(rows,1): 
            f.write(f"{rank} {-math.log(lam):.16e} {lam:.16e} {q} {iq}\n")
    return rows



def clusters(x,tol):
    """Group neighbouring entanglement levels within an absolute tolerance.

    This is a simple practical definition of a quasi-degenerate multiplet.
    Change --tol when the expected level splitting or numerical accuracy changes.
    """
    c=[]
    for v in x:
        if not c or v-c[-1][-1]>tol: 
            c.append([v])
        else: c[-1].append(v)
    return c

def statistics(rows):
    """Return symmetry-sector and Schmidt-spectrum diagnostics.

    P_q is the normalized weight of sector q.  Sn is the entropy of the
    sector-weight distribution, while Sc is the entropy internal to sectors.
    Their sum is the total von Neumann entropy.  D_eff=exp(Stotal) is the
    exponential Schmidt rank; IPR=sum(lambda**2) is the purity, and N_part is
    its inverse.  The spectrum is renormalized to the retained weight, which
    is useful when DMRG discarded a small tail of eigenvalues.
    """
    byq={}
    for lam,q,iq in rows: 
        byq[q]=byq.get(q,0.0)+lam
    norm=sum(byq.values()); 
    weights={q:w/norm for q,w in byq.items()}
    sn=-sum(w*math.log(w) for w in weights.values() if w>0); sc=0.0
    for lam,q,iq in rows:
        p=(lam/norm)/weights[q]
        if p>0: sc-=weights[q]*p*math.log(p)
    purity=sum((lam/norm)**2 for lam,q,iq in rows)
    return weights,sn,sc,sn+sc,math.exp(sn+sc),purity,1/purity

def main():
    """Run extraction, write tables, and make the three diagnostic figures."""
    ap=argparse.ArgumentParser(); 
    # Number of low-lying xi levels shown in the spectrum plot.
    # --tol is the criterion used to identify quasi-degenerate multiplets.
    ap.add_argument("--root",type=Path,default="./")
    ap.add_argument("--nlevels",type=int,default=20) 
    ap.add_argument("--tol",type=float,default=2e-10)
    #
    a=ap.parse_args(); 
    #
    out=a.root/"entanglement_spectrum"; 
    out.mkdir(exist_ok=True)
    #
    # The glob is deliberately independent of L, so the same script works for
    # L28, L40, ... and for every M* directory with the same data layout.
    files=sorted(a.root.glob("U*/LambdaQ_*_L*_*.dmrg"))
    if not files: 
        raise SystemExit("No LambdaQ files found")
    agg={"left":[],"right":[]}; 
    gaps=[]; 
    multi=[]; 
    stats={"left":[],"right":[]}; 
    weights={"left":{},"right":{}}
    #
    # Read both RDM sides.  The left side is also used for the gap plots;
    # the right side is retained for consistency checks and diagnostics.
    for p in files:
        side="left" if "_left_" in p.name else "right"; 
        u=float(p.parent.name[1:]); 
        rows=extract(p,side,out)
        #
        w,sn,sc,st,deff,ipr,npart=statistics(rows); 
        #
        stats[side].append((u,sn,sc,st,deff,ipr,npart)); 
        weights[side][u]=w
        #
        agg[side].extend((u,i,-math.log(x),x,q,iq) for i,(x,q,iq) in enumerate(rows,1))
        if side=="left" and len(rows)>=2:
            xi=[-math.log(x[0]) for x in rows]; 
            gaps.append((u,xi[0],xi[1],xi[1]-xi[0],rows[0][0],rows[1][0]))
            c=clusters(xi[:a.nlevels],a.tol); 
            multi.append((u,c[0][0],c[1][0],c[1][0]-c[0][0],len(c[0]),len(c[1])))
    # Raw, gnuplot-friendly aggregate spectra.
    for side,data in agg.items():
        with (out/f"spectrum_{side}_all.dat").open("w") as f:
            f.write("# U rank xi lambda q index_in_sector\n"); 
            [f.write(" ".join(map(str,r))+"\n") for r in sorted(data)]
    # Gaps between individual levels and between quasi-degenerate multiplets.
    with (out/"entanglement_gap.dat").open("w") as f:
        f.write("# U xi1 xi2 gap lambda1 lambda2\n"); 
        [f.write(" ".join(f"{x:.12e}" for x in r)+"\n") for r in sorted(gaps)]
    # Gaps between quasi-degenerate multiplets.       
    with (out/"multiplet_gaps.dat").open("w") as f:
        f.write("# U xi_m1 xi_m2 gap_multiplets size_m1 size_m2\n"); 
        [f.write(" ".join(f"{x:.12e}" for x in r)+"\n") for r in sorted(multi)]
    # Sector weights and scalar Schmidt diagnostics, one row per U.
    for side in ("left","right"):
        with (out/f"schmidt_statistics_{side}.dat").open("w") as f:
            f.write("# U Sn Sc Stotal D_eff IPR N_participation\n")
            [f.write(" ".join(f"{x:.12e}" for x in r)+"\n") for r in sorted(stats[side])]
        qvals=sorted({q for w in weights[side].values() for q in w})
        with (out/f"sector_weights_{side}.dat").open("w") as f:
            f.write("# U "+" ".join(f"P_q{q}" for q in qvals)+"\n")
            for u in sorted(weights[side]): f.write(f"{u:.12e} "+" ".join(f"{weights[side][u].get(q,0):.12e}" for q in qvals)+"\n")
    #
    #
    # Figure 1: low-lying entanglement spectrum and multiplet gap.
    fig,(ax0,ax1)=plt.subplots(2,1,figsize=(10,8),sharex=True)
    #
    for p in sorted(out.glob("spectrum_left_U*.dat")):
        u=float(p.stem.split("_U")[1]); 
        d=np.loadtxt(p,comments="#"); 
        d=np.atleast_2d(d)[:a.nlevels]
        ax0.scatter([u]*len(d),d[:,1],c=d[:,3],cmap="tab10",s=13,vmin=0,vmax=10)
    #
    #       
    ax0.set_ylim(0,12); 
    ax0.set_ylabel(r"$\xi=-\log\lambda$"); 
    ax0.set_title(f"BHZ entanglement spectrum: {a.root.name}"); 
    ax0.grid(alpha=.25)
    #
    #
    g=np.array(multi); 
    ax1.plot(g[:,0],g[:,3],"o-",label="gap between multiplets"); 
    ax1.plot(g[:,0],g[:,3]*0+np.array([r[3] for r in gaps[:len(g)]]),"s--",label="first-level gap")
    ax1.set_xlabel("U"); 
    ax1.set_ylabel(r"$\Delta\xi$"); 
    ax1.set_ylim(bottom=0); 
    ax1.grid(alpha=.25); 
    ax1.legend(frameon=False)
    fig.tight_layout(); 
    fig.savefig(out/"entanglement_evolution.png",dpi=300)
    #
    # Figure 2: entropy decomposition and effective Schmidt complexity.
    fig,axs=plt.subplots(2,2,figsize=(11,8),sharex=True)
    for side,color in (("left","#245b9c"),("right","#b33b3b")):
        s=np.array(sorted(stats[side])); 
        axs[0,0].plot(s[:,0],s[:,1],"o-",color=color,label=side+" $S_n$"); 
        axs[0,0].plot(s[:,0],s[:,2],"--",color=color,label=side+" $S_c$")
        axs[0,1].plot(s[:,0],s[:,4],"o-",color=color,label=side+" $D_{eff}$"); 
        axs[1,0].plot(s[:,0],s[:,5],"o-",color=color,label=side+" IPR"); 
        axs[1,1].plot(s[:,0],s[:,6],"o-",color=color,label=side+" $N_{part}$")
    for ax in axs.flat: 
        ax.grid(alpha=.25); 
    ax.legend(frameon=False,fontsize=8)
    axs[0,0].set_ylabel("entropy"); 
    axs[0,1].set_ylabel(r"$e^{S_{tot}}$"); 
    axs[1,0].set_ylabel(r"$\sum_i\lambda_i^2$"); 
    axs[1,1].set_ylabel(r"$1/\sum_i\lambda_i^2$"); 
    axs[1,0].set_xlabel("U"); 
    axs[1,1].set_xlabel("U")
    fig.tight_layout(); 
    fig.savefig(out/"schmidt_statistics.png",dpi=300)


    # Figure 3: heatmaps of P_q(U), separately for left and right RDMs.
    fig,axs=plt.subplots(2,1,figsize=(10,7),sharex=True)
    for ax,side in zip(axs,("left","right")):
        qvals=sorted({q for w in weights[side].values() for q in w}); 
        us=sorted(weights[side]); 
        z=np.array([[weights[side][u].get(q,0) for u in us] for q in qvals])
        im=ax.pcolormesh(us,qvals,z,shading="nearest",cmap="viridis",vmin=0,vmax=max(1e-12,z.max())); 
        ax.set_ylim(1,10)
        ax.set_ylabel("q"); 
        ax.set_title(f"Sector weights $P_q(U)$ ({side})"); 
        fig.colorbar(im,ax=ax,label="$P_q$")

    axs[-1].set_xlabel("U"); 
    fig.tight_layout(); 
    fig.savefig(out/"sector_weights_heatmap.png",dpi=300)
    print(f"Wrote analysis to {out}")




    
if __name__=="__main__": main()
