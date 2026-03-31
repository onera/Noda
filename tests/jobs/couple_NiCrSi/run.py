import numpy as np

from noda import simu

simu1 = simu.NewSimulation(file='couple_NiCrSi.toml')
simu1.run()

#%% Plot

fig, ax = simu1.plot()

#%% Validation

res = simu1.result()
x_simu = np.array(list(res.x.values()))

ref = np.genfromtxt('couple_NiCrSi-ref.txt', skip_header=1, delimiter=',')
x_ref = ref.T[1:]

assert np.allclose(x_simu, x_ref, atol=1e-6)
