import numpy as np

from noda import simu

foo = simu.NewSimulation('couple_NiCrSi')
foo.run()

#%% Plot

fig, ax = foo.plot()

#%% Validation

bar = foo.result()
x_simu = np.array(list(bar.x.values()))

ref = np.genfromtxt('couple_NiCrSi-ref.txt', skip_header=1, delimiter=',')
x_ref = ref.T[1:]

assert np.allclose(x_simu, x_ref, atol=1e-6)
