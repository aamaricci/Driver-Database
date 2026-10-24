"""Analytic atomic spectra and optional regression against an older solver."""
from pathlib import Path
import argparse
import subprocess
import tempfile
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--exe', type=Path, default=Path(__file__).with_name('atomic_bhz'))
parser.add_argument('--reference', type=Path)
args = parser.parse_args()
exe = args.exe.resolve()
reference = args.reference.resolve() if args.reference else None
worst = 0.0
with tempfile.TemporaryDirectory(prefix='check_atomic_') as work:
    root = Path(work)
    for n in range(5):
        for full in (0, 1):
            for hf in (False, True):
                for ratio in (-1, .25):
                    config = root / 'inputATOM.conf'
                    config.write_text(
                        f'NELEC={n}\nNU=5\nNM=9\nUMIN=0\nUMAX=8\n'
                        f'MMIN=-2\nMMAX=2\nJH=.5\nJRATIO={ratio}\n'
                        f'JX_FACTOR={full}\nJP_FACTOR={full}\nHFMODE={"T" if hf else "F"}\n')
                    def run(binary, folder):
                        folder.mkdir(exist_ok=True)
                        (folder / 'inputATOM.conf').write_text(config.read_text())
                        with (folder / 'run.log').open('w') as log:
                            subprocess.run([str(binary), 'FINPUT=inputATOM.conf'],
                                           cwd=folder, stdout=log, stderr=log, check=True)
                        return tuple(np.loadtxt(folder / name) for name in
                                     ('atomic_ground.dat', 'atomic_spectrum.dat'))
                    ground, spectrum = run(exe, root / 'new')
                    assert np.allclose(ground[:, -2:].sum(axis=1), n, atol=1e-12), (n,full,hf,ratio,ground[:2])
                    for row in ground:
                        u, m, j = row[:3]
                        levels = {
                            0: [0],
                            1: [-m, -m, m, m],
                            2: [u-3*j, u-3*j, u-2*j-full*j, u-2*j+full*j,
                                u-np.hypot(2*m, full*j), u+np.hypot(2*m, full*j)],
                            3: [3*u-5*j-m]*2 + [3*u-5*j+m]*2,
                            4: [6*u-10*j],
                        }[n]
                        expected = np.sort(levels) - (n*(3*u-5*j)/2 if hf else 0)
                        actual = spectrum[(spectrum[:, 0] == u) &
                                          (spectrum[:, 1] == m), -1]
                        assert np.allclose(actual, expected, atol=1e-12, rtol=0)
                        assert abs(row[3] - expected[0]) < 1e-12
                        g = np.count_nonzero(abs(expected-expected[0]) <
                                             1e-10*max(1, max(abs(expected))))
                        assert int(row[4]) == g
                        if n == 2 and g == (3 if full else 2):
                            if u-3*j < u-np.hypot(2*m, full*j)-1e-10:
                                assert abs(row[6]-2) < 1e-12
                                assert abs(row[8]) < 1e-12
                        if n == 2 and g == 1:
                            assert abs(row[6]) < 1e-12
                            assert abs(row[9]-1) < 1e-12
                            assert abs(row[8]+2*m/np.hypot(2*m, full*j)) < 1e-12
                    if reference:
                        old = run(reference, root / 'old')
                        for new_data, old_data in zip((ground, spectrum), old):
                            delta = np.max(abs(new_data-old_data))
                            worst = max(worst, delta)
                            assert np.allclose(new_data, old_data, atol=1e-12, rtol=0)
print('PASS: analytic spectra, degeneracies, spin, orbital polarization; '
      'N=0..4, density-density/full Kanamori, HFMODE on/off, fixed J and J/U.')
if reference:
    print(f'Maximum difference against reference: {worst:.3e}')
