import numpy as np
import matplotlib.pyplot as plt

from noda import simu

def analytical_solution(r, R, D, ts, x_bulk, x_surf, N=100):
    f = 1 + 2*sum(((-1)**n)
                  * np.divide(R*np.sin(n*np.pi*r/R), np.pi*n*r,
                              out=np.ones(r.size), where=(r!=0))
                  * np.exp(-D*n**2*np.pi**2*ts/R**2)
                  for n in range(1, N))
    x = x_bulk + (x_surf - x_bulk)*f
    return x

s = simu.NewSimulation('sphere_AB')
s.run()

#%% Comparison with analytical solution

x_bulk = s.x_init['B'][0]
x_surf = s.BC['c_right'](0).x.mid[0]
R = s.zmax
r = s.z_init

x_any = np.array([[0.5]])
D = s.DT_fun(x_any)[0]

th_targets = [1, 5, 10]
colors = ['salmon', 'gold', 'paleturquoise']

fig, ax = plt.subplots()
for th_target, c in zip(th_targets, colors):
    th = s.saved_times[np.argmin(abs(th_target - s.saved_times))]
    res = s.result(th=th)
    x_ana = analytical_solution(r, R, D, th*3600, x_bulk, x_surf)
    ax.plot(r*1e6, x_ana, 'k')
    ax.plot(r*1e6, res.x['B'], '--', c=c, label=f"{th:.0f}")
ax.set_xlabel(r'radius ($\mu$m)')
ax.set_ylabel('$x_B$')
ax.legend(title='time (h)')
ax.grid(visible=True)
fig.set_dpi(200)

# plt.savefig('sphere_AB.png', bbox_inches="tight", dpi=100)

#%% Validation

assert np.allclose(res.x['B'], x_ana, atol=2e-3)
