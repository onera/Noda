import numpy as np

from noda import simu

foo = simu.AlloySystem('NiCrSi')

# The independent constituents are Cr and Si, in alphabetical order:
assert foo.inds == ['Cr', 'Si']

# Make atom fraction array
# The required shape is (n_inds, n_points), where n_inds is the number of
# independent constituents, and n_points the number of compositions of
# interest. Here n_inds = 2. The independent constituents must be entered in
# the order used by foo.inds (alphabetical).
x_Cr = np.array([0.1, 0.2, 0.3])
x_Si = np.array([0.1, 0.1, 0.01])
x = np.vstack((x_Cr, x_Si))

# Compute Gibbs free energy, shape = (n_points,)
G = foo.G_fun(x)

# Compute chemical potentials, shape = (n ind. constituents, n_points)
MU = foo.MU_fun(x)

# Compute tracer diffusion coefficients, shape = (n ind. constituents, n_points)
DT = foo.DT_fun(x)


#%% Validation

ref = np.genfromtxt('NiCrSi-ref.txt', skip_header=1, delimiter=',').T
x_ref = ref[:2]
G_ref = ref[2:3]
MU_ref = ref[3:6]
DT_ref = ref[6:9]

assert np.allclose(x, x_ref)
assert np.allclose(G, G_ref)
assert np.allclose(MU, MU_ref)
assert np.allclose(DT, DT_ref)
