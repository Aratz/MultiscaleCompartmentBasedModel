# WMM
import sys
import WMM
import time
import reactionrates as rr

def adjust_parameters(parameters):
    parameters = {**parameters}
    sigma = 0.01
    V_nucleus = 4*np.pi*parameters['nucleus_radius']**3/3
    k_a_wmmodel, k_d_wmmodel = rr.well_mixed_rates(
        parameters['k_a'],
        parameters['k_d'],
        sigma, 2*parameters['diffusion'],
        V_nucleus)

    parameters['k_a'] = k_a_wmmodel
    parameters['k_d'] = k_d_wmmodel

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

init_model = WMM.Cell(adjust_parameters(base_parameters))

# Generate initial state
initial_state = [smoldyn.Mol("Gf", 0.0, 0.0, 0.0)]
_, initial_state = init_model.run(initial_state, time_stop=INIT_STOP, n_trajectories=1, n_points=2)

wmm_parameters = {**base_parameters}
wmm_parameters['mu'] *= 2**chi
wmm_parameters['kappa'] *= 2**chi
wmm_parameters['gamma'] *= 2**chi
wmm_parameters['diffusion'] *= 2**D
wmm_parameters['k_d'] *= 2**k_d

model = WMM.Cell(adjust_parameters(wmm_parameters))

start = time.time()
history = model.run(initial_state, time_stop=TSTOP, n_trajectories=NTRAJ, n_points=500)
end = time.time()

with open('data/wmm_data.json', 'w') as f:
json.dump((wmm_parameters, history), f)

with open('data/wmm_time.json', 'w') as f:
json.dump((wmm_parameters, end - start), f)
