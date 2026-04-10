import pytest
import numpy as np

from noda.simu import NewSimulation
from noda.utils import UserInputError

config = {'databases': {'thermo': "AB_thermo_ideal",
                        'mobility': "XXX"},
          'temperature': {'TC': 1200},
        }

# Missing 'system' table
with pytest.raises(UserInputError):
    NewSimulation(config=config, log=False)

config['system'] = {'components': 'A, B, C', 'phases': 'fcc'}

# Element C not found in thermodynamic database
with pytest.raises(UserInputError):
    NewSimulation(config=config, log=False)

config['system']['components'] = 'A, B'

# Mobility database not found in 'user_data.toml'
with pytest.raises(UserInputError):
    NewSimulation(config=config, log=False)

config['databases']['mobility'] = 'AB_mob_ideal'

VA = 5e-6
G_Va_A = [2e5, 60]
G_Va_B = [1e5, 50]
config['databases']['partial_molar_volume'] = {'A': VA}
config['databases']['vacancy_formation_energy'] = {'A': G_Va_A, 'B': G_Va_B}
x = np.array([[0.1, 0.2, 0.3]])
yVa = np.array([0.00022117, 0.0004437 , 0.0008901 ])

s = NewSimulation(config=config, log=False)
assert s.V_partial['A'] == VA
assert np.allclose(s.thermo.yVa_fun(x), yVa)

config['initial_conditions'] = {'atom_fraction': {'shape': 'some_shape',
                                                  'step_fraction': 0.5,
                                                  'right': {'B': 0.3}
                                                  }
                                }

# Missing 'space' table in configuration
with pytest.raises(UserInputError):
    NewSimulation(config=config, log=False)

config['space'] = {'zmax': 6e-4, 'nz': 100}

# Initial conditions : invalid atom fraction shape
with pytest.raises(UserInputError):
    NewSimulation(config=config, log=False)

config['initial_conditions']['atom_fraction']['shape'] = 'step'

# Initial conditions : missing 'left' entry
with pytest.raises(UserInputError):
    NewSimulation(config=config, log=False)

config['initial_conditions']['atom_fraction']['left'] = {'B': 0.1, 'C': 0.5}

# Initial conditions : extra component on left side
with pytest.raises(UserInputError):
    NewSimulation(config=config, log=False)

config['initial_conditions']['atom_fraction']['left'].pop('C')

# Missing time table
with pytest.raises(UserInputError):
    s = NewSimulation(config=config)
    s.run()

config['time'] = {}

# Time : require either th or ts
with pytest.raises(UserInputError):
    NewSimulation(config=config, log=False)

config['time'] = {'th': 10, 'ts': 2}

# Time : require either th or ts
with pytest.raises(UserInputError):
    NewSimulation(config=config, log=False)

config['time'] = {'th': 10}
k_dislo = 0
config['options'] = {'k_dislo' : k_dislo, 'rho_dislo' : 0}

# Cannot specify both k_dislo and rho_dislo.
with pytest.raises(UserInputError):
    NewSimulation(config=config, log=False)

config['options'].pop('rho_dislo')
s = NewSimulation(config=config, log=False)
assert s.lattice.k_dislo == k_dislo
