"""Compare DMRG against full-basis spin-chain ED in OBC and PBC."""
from pathlib import Path
import subprocess
import sys
import numpy as np

exe = str(Path(sys.argv[1]).resolve())
root = Path(sys.argv[2]).resolve()
root.mkdir(parents=True, exist_ok=True)
basis = [s for s in range(64) if s.bit_count() == 3]
lookup = {s: i for i, s in enumerate(basis)}

def gs(periodic):
    h = np.zeros((len(basis), len(basis)))
    for col, s in enumerate(basis):
        for p in range(6 if periodic else 5):
            n1, n2 = (s >> p) & 1, (s >> ((p+1)%6)) & 1
            h[col, col] += (n1-.5)*(n2-.5)
            if n1 != n2:
                h[lookup[s ^ (1 << p) ^ (1 << ((p+1)%6))], col] += .5
    return np.linalg.eigh(h)[1][:, 0]

for periodic in (False, True):
    work = root / ('pbc' if periodic else 'obc')
    work.mkdir(exist_ok=True)
    (work / 'DMRG.conf').write_text(f'''DMRGTYPE=i
    PBCDMRG={'T' if periodic else 'F'}
    LDMRG=3
    MDMRG=256
    EDMRG=0
    QNDIM=1
    QNTYPE=local
    DMRG_QN=0
    NORB=1
    LANC_NEIGEN=1
    LANC_DIM_THRESHOLD=4096
    SPARSE_H=F
    SAVE_BLOCK=T
    SAVE_UMAT=T
    BLOCK_UMAT_CACHE=T
    SAVE_MEASURE_STATE=T
    RESTART=F
    ''')
    command = [exe]
    if len(sys.argv) > 3:
        command = [sys.argv[3], '-n', '2', exe]
    result = subprocess.run(command, cwd=work, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=120)
    (work / 'run.log').write_text(result.stdout)
    assert result.returncode == 0, result.stdout[-6000:]
    state = gs(periodic)
    polarization = np.array([sum((p+1)*((s >> p)&1) for p in range(6)) for s in basis])
    expected = [np.sum(abs(state)**2*np.exp(1j*2*np.pi*scale*polarization/6)) for scale in (1, 2)]
    data = list(map(float, (work / 'polarization.out').read_text().split()))
    actual = [complex(data[0], data[1]), complex(data[2], data[3])]
    assert np.max(np.abs(np.array(actual)-expected)) < 2e-9, (periodic, actual, expected)
    print(f"{'PBC' if periodic else 'OBC'}: DFT and complex marker passed against ED")
