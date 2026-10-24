from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

folder = Path(__file__).resolve().parent
print(folder)
data = np.loadtxt(folder / 'atomic_ground.dat')
u = np.unique(data[:, 0])
m = np.unique(data[:, 1])
s2 = data[:, 6].reshape(len(u), len(m)).T
uu, mm = np.meshgrid(u, m)
boundary = np.isclose(mm, 3*uu/8, atol=1e-10, rtol=0)
assert np.allclose(s2[(mm < 3*uu/8) & ~boundary], 2)
assert np.allclose(s2[(mm > 3*uu/8) & ~boundary], 0)
phase = (s2 > 1).astype(int)

fig, ax = plt.subplots(figsize=(8.1, 6.1), constrained_layout=True)

ax.pcolormesh(u, m, phase, shading='nearest', rasterized=True,
              cmap=ListedColormap(['#eac687', '#80b3cc']), vmin=0, vmax=1)
ax.plot(u, 3*u/8, color='#263643', lw=2,
        label=r'$M_c=\frac{3}{2}J=\frac{3}{8}U$')
ax.text(4.0, 3.25, r'$n_1=0,\ n_2=2,\quad S^2=0$, $g=1$',
        ha='center', va='center', fontsize=12, linespacing=1.5)
ax.text(5.5, 0.65,
        r'$n_1=n_2=1,\quad S_z=\pm1$, $S^2=2$, $g=2$',
        ha='center', va='center', fontsize=12, linespacing=1.5)
ax.set(xlim=(0, 8), ylim=(0, 4), xlabel=r'$U$', ylabel=r'$M$')
ax.set_title(r'BHZ $N=2$, $t=\lambda=0,\quad J/U=0.25,\quad J_x=J_p=0$',
             fontsize=14, pad=13)
ax.tick_params(labelsize=11)
ax.xaxis.label.set_size(14)
ax.yaxis.label.set_size(14)
ax.legend(loc='upper left', framealpha=.95, fontsize=12)
fig.get_layout_engine().set(rect=(0, .055, 1, 1))
# fig.text(.5, .018, 'Sulla frontiera: degenerazione 3 per U > 0. '
#          'Colori da diagonalizzazione esatta su griglia 321 × 321.',
#          ha='center', fontsize=9, color='#333333')
out = folder / 'phase_diagram_U_M_J025_density_density.png'
fig.savefig(out, dpi=220, facecolor='white', bbox_inches='tight', pad_inches=.16)
print(out)
