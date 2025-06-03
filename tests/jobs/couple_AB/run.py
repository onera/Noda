import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt

from noda import simu

s = simu.NewSimulation('couple_AB')
s.run()

#%% Comparison with analytical solution

x_left = s.x_init['B'][0]
x_right = s.x_init['B'][-1]

x_any = np.array([[0.5]])
D = s.DT_fun(x_any)[0]

r = s.results[-1]
zmid = (r.z[-1] + r.z[0])/2
x_ana = x_right + (x_left - x_right)*0.5*erfc((r.z - zmid)/np.sqrt(4*D*s.ts))

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1,
                               figsize=(5, 5),
                               gridspec_kw={'height_ratios': [3, 1]})
ax1.plot(r.z*1e6, x_ana, 'k', label='analytical solution')
ax1.plot(r.z*1e6, r.x['B'], '--', c='tab:orange', label='Noda')
ax1.set_ylabel('$x_B$')
ax1.set_xticklabels([])
ax1.legend()
ax1.grid(visible=True)

ax2.plot(r.z*1e6, r.x['B'] - x_ana, c='C1')
ax2.set_xlabel(r'$z$ ($\mu$m)')
ax2.set_ylabel(r'$x_B^\mathrm{Noda} - x_B^\mathrm{ana}$')
ax2.grid(visible=True)

# plt.savefig('couple_AB.png', bbox_inches="tight", dpi=100)

#%% Validation

assert np.allclose(r.x['B'], x_ana, atol=2e-4)
