import numpy as np
import matplotlib.pyplot as plt

from noda import simu
from noda.constants import R

s = simu.AlloySystem('NiCr')

#%% Prepare data

comps = [k for k in s.comps if k != 'Va']

x0 = np.linspace(1e-9, 0.4, num=200)
x = x0[np.newaxis]

DT = s.DT_fun(x)
MU = s.MU_fun(x)
phi = -(1 - x0)/(R*s.TK)*np.gradient(MU[1], x0)
DI = DT*phi

#%% Plot

colors = ['tab:blue', 'tab:orange']

fig, ax = plt.subplots(figsize=(5,5/1.5))
for i, k in enumerate(comps):
    ax.plot(x0, np.log10(DT[i]), c=colors[i], label=k)
    ax.plot(x0, np.log10(DI[i]), '--', c=colors[i])

ax.set_xlabel(r'$x_\mathrm{Cr}$')
ax.set_ylabel(r'$\log_{10}\ D_k$ (m$^2$/s)')
ax.set_ylim(-14.2, -13.2)
fs = 14
xpos = 0.25
ax.annotate(r'$\bar{D}_\mathrm{Ni}$', (xpos, -13.92), fontsize=fs)
ax.annotate(r'$D^*_\mathrm{Ni}$', (xpos, -14.16), fontsize=fs)
ax.annotate(r'$\bar{D}_\mathrm{Cr}$', (xpos, -13.48), fontsize=fs)
ax.annotate(r'$D^*_\mathrm{Cr}$', (xpos, -13.72), fontsize=fs)
# plt.savefig('diffusivities_NiCr.png', bbox_inches="tight", dpi=100)
