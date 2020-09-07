import sys
import CBM
import time
import json
import numpy as np
import reactionrates as rr

def adjust_parameters(parameters):
    parameters = {**parameters}
    kentry, kexit = rr.entry_exit_rates(
        parameters['nucleus_radius'],
        parameters['cell_radius'],
        parameters['diffusion']) 

    sigma = 0.01
    V_nucleus = 4*np.pi*parameters['nucleus_radius']**3/3
    k_a_cbmodel, k_d_cbmodel = rr.well_mixed_rates(
        parameters['k_a'],
        parameters['k_d'],
        sigma, 2*parameters['diffusion'],
        V_nucleus)

    parameters['k_a'] = k_a_cbmodel
    parameters['k_d'] = k_d_cbmodel
    parameters['k_cn'] = kentry
    parameters['k_nc'] = kexit

    return parameters

base_parameters = {
    'k_a': 0.00166,
    'k_d': 0.1,
    'mu': 3.0,
    'kappa': 1.0,
    'gamma': 0.04,
    'diffusion': 0.6,
    'time_step': 0.32,
    'cell_radius': 6.0,
    'nucleus_radius': 2.5,
}

INIT_STOP = 200
TSTOP = 1200 # The 200 first minutes will be discarded ad burn-in
NTRAJ = 64
NPOINTS = 120
THRES = 0.1

D, chi, k_d = [float(v) for v in sys.argv[1:4]]

init_model = CBM.Cell(adjust_parameters(base_parameters))

# Generate initial state
initial_state = {'Gf': 1}
_, initial_state = init_model.run(initial_state, time_stop=INIT_STOP, n_trajectories=1, n_points=2)

cbm_parameters = {**base_parameters}
cbm_parameters['mu'] *= 2**chi
cbm_parameters['kappa'] *= 2**chi
cbm_parameters['gamma'] *= 2**chi
cbm_parameters['diffusion'] *= 2**D
cbm_parameters['k_d'] *= 2**k_d

model = CBM.Cell(adjust_parameters(cbm_parameters))

start = time.time()
history = model.run(initial_state, time_stop=TSTOP, n_trajectories=NTRAJ, n_points=NPOINTS)
end = time.time()

with open(f'data/CBM_data({D},{chi},{k_d}).json', 'w') as f:
    json.dump((cbm_parameters, history), f)

with open(f'data/CBM_time({D},{chi},{k_d}).json', 'w') as f:
    json.dump((cbm_parameters, end - start), f)
