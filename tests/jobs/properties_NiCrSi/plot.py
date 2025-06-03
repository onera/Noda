import numpy as np
import matplotlib.pyplot as plt

from noda import simu

s = simu.AlloySystem('NiCrSi')

#%% Gibbs free energy

fig, ax = plt.subplots()
for val in [0.01, 0.05, 0.1]:
    x0 = np.linspace(1e-3, 0.4)
    x1 = np.ones(x0.size)*val
    x = np.vstack((x0, x1))
    G = s.G_fun(x)
    ax.plot(x0, G/1000, label=f'{val:.2f}')
ax.set_xlabel(r'$x_\mathrm{Cr}$')
ax.set_ylabel('$G$ (kJ/mol)')
ax.legend(title=r'$x_\mathrm{Si}$', loc='upper left')

#%% Chemical potentials

fig, axes = plt.subplots(3, 2, figsize=(6, 8), tight_layout=True)

for i in range(3):
    
    A = s.comps[i + 1]
    ax1 = axes[i, 0]
    ax2 = axes[i, 1]

    targets = [0.01, 0.05, 0.1]
    x0 = np.linspace(1e-3, 0.4, num=200)
    
    for j, val in enumerate(targets):
        x1 = np.ones(x0.size)*val
        x = np.vstack((x0, x1))
        MU = s.MU_fun(x)[i]/1000
        ax1.plot(x0, MU, label=f'{val:.2f}')
    
    ax1.set_ylabel(rf'$\mu_\mathrm{{{A}}}$ (kJ/mol)')
    ax1.legend(title=r'$x_\mathrm{Si}$', loc='upper left')
    
    targets = [0.1, 0.2, 0.3]
    x1 = np.linspace(1e-3, 0.12, num=200)
    
    for j, val in enumerate(targets):
        x0 = np.ones(x1.size)*val
        x = np.vstack((x0, x1))
        MU = s.MU_fun(x)[i]/1000
        ax2.plot(x1, MU, label=f'{val:.2f}')
    
    ax2.legend(title=r'$x_\mathrm{Cr}$', loc='upper left')

ax1.set_xlabel(r'$x_\mathrm{Cr}$')
ax2.set_xlabel(r'$x_\mathrm{Si}$')

# plt.savefig('chemical_potentials_NiCrSi.png', bbox_inches="tight", dpi=100)

#%% Tracer diffusion coefficients

fig, axes = plt.subplots(3, 2, figsize=(6, 8), tight_layout=True)

for i in range(3):
    
    A = s.comps[i + 1]
    ax1 = axes[i, 0]
    ax2 = axes[i, 1]

    targets = [0.01, 0.05, 0.1]
    x0 = np.linspace(1e-3, 0.4, num=200)
    
    for j, val in enumerate(targets):
        x1 = np.ones(x0.size)*val
        x = np.vstack((x0, x1))
        DT = s.DT_fun(x)[i]
        ax1.plot(x0, np.log10(DT), label=f'{val:.2f}')
    
    ax1.set_ylabel(rf'$\log_{{{10}}}\ D_\mathrm{{{A}}}^*$ (m$^2$/s)')
    ax1.legend(title=r'$x_\mathrm{Si}$', loc='upper left')
    
    targets = [0.1, 0.2, 0.3]
    x1 = np.linspace(1e-3, 0.12, num=200)
    
    for j, val in enumerate(targets):
        x0 = np.ones(x1.size)*val
        x = np.vstack((x0, x1))
        DT = s.DT_fun(x)[i]
        ax2.plot(x1, np.log10(DT), label=f'{val:.2f}')
    
    ax2.legend(title=r'$x_\mathrm{Cr}$', loc='upper left')

ax1.set_xlabel(r'$x_\mathrm{Cr}$')
ax2.set_xlabel(r'$x_\mathrm{Si}$')
