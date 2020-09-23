"""
Usage: python3 run_inference_CBM.py <data.json> <output.db>

"""

import sys

import json
import pyabc
import numpy as np
import scipy.stats as st

import CBM
import WMM
import reactionrates as rr

if len(sys.argv) <= 1:
    raise Exception("Input file missing")
if len(sys.argv) <= 2:
    raise Exception("Output file missing")

EPSILON = 10
MAXPOP = 10


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

def adjust_parameters(parameters, model):
    print(parameters)
    parameters = {**parameters}

    sigma = 0.01
    V_nucleus = 4*np.pi*parameters['nucleus_radius']**3/3
    k_a_cbmodel, k_d_cbmodel = rr.well_mixed_rates(
        parameters['k_a'],
        parameters['k_d'],
        sigma, 2*parameters['diffusion'],
        V_nucleus)
    parameters['k_a'] = k_a_cbmodel
    parameters['k_d'] = k_d_cbmodel

    if model == 'CBM':
        kentry, kexit = rr.entry_exit_rates(
            parameters['nucleus_radius'],
            parameters['cell_radius'],
            parameters['diffusion']) 
        parameters['k_cn'] = kentry
        parameters['k_nc'] = kexit

    return parameters

INIT_STOP = 200
TSTOP = 1200 # The 200 first minutes will be discarded ad burn-in
NTRAJ = 64
NPOINTS = 120
THRES = 0.1
BURNIN = 20

def neg_feedback_WMM(params):
    wmm_parameters = {**base_parameters}
    wmm_parameters['mu'] *= 2**params['chi']
    wmm_parameters['kappa'] *= 2**params['chi']
    wmm_parameters['gamma'] *= 2**params['chi']
    wmm_parameters['diffusion'] *= 2**params['diffusion']
    #wmm_parameters['k_d'] *= 2**params['k_d']

    init_model = WMM.Cell(adjust_parameters(base_parameters, 'WMM'))

    # Generate initial state
    initial_state = {'Gf': 1}
    _, initial_state = init_model.run(initial_state, time_stop=INIT_STOP, n_trajectories=1, n_points=2)

    model = WMM.Cell(adjust_parameters(wmm_parameters, 'WMM'))

    results = model.run({'Gf':1}, time_stop=TSTOP, n_trajectories=NTRAJ, n_points=NPOINTS)

    return {
            "P": np.array(results['P'])[:, BURNIN:].flatten().tolist(),
            "RNA": np.array(results['RNA'])[:, BURNIN:].flatten().tolist(),
        }

def neg_feedback_CBM(params):
    cbm_parameters = {**base_parameters}
    cbm_parameters['mu'] *= 2**params['chi']
    cbm_parameters['kappa'] *= 2**params['chi']
    cbm_parameters['gamma'] *= 2**params['chi']
    cbm_parameters['diffusion'] *= 2**params['diffusion']
    #cbm_parameters['k_d'] *= 2**params['k_d']

    init_model = CBM.Cell(adjust_parameters(base_parameters, 'CBM'))

    # Generate initial state
    initial_state = {'Gf': 1}
    _, initial_state = init_model.run(initial_state, time_stop=INIT_STOP, n_trajectories=1, n_points=2)

    model = CBM.Cell(adjust_parameters(cbm_parameters, 'CBM'))

    results = model.run({'Gf':1}, time_stop=TSTOP, n_trajectories=NTRAJ, n_points=NPOINTS)

    return {
            "P": (np.array(results['Pnuc']) + np.array(results['Pcyt']))[:, BURNIN:].flatten().tolist(),
            "RNA": (np.array(results['RNAnuc']) + np.array(results['RNAcyt']))[:, BURNIN:].flatten().tolist(),
        }

models = [neg_feedback_WMM, neg_feedback_CBM]

bounds = {
    "chi": (-2, 4),
    "diffusion": (-8, 4),
    #"k_d": (-3, 3),
}

priors = [
    pyabc.Distribution(
        **{param: pyabc.RV("uniform", p_min, p_max - p_min)
        #**{param: pyabc.RV("loguniform", p_min, p_max)
          for param, (p_min, p_max) in bounds.items()}
    )
    for _ in models
]

summary_statistics_raw = {
        'mean': lambda x: np.mean(x),
        'std': lambda x: np.std(x),
        'min': lambda x: np.min(x),
        'max': lambda x: np.max(x),
        }

stat_mean = {'meanP': 1042.154285695293, 'meanRNA': 61.844218439788314, 'stdP': 269.7571445476583, 'stdRNA': 15.925807595648152, 'minP': 242.33203125, 'minRNA': 11.462890625, 'maxP': 1731.669921875, 'maxRNA': 104.96223958333333}
stat_std = {'meanP': 526.6189072354466, 'meanRNA': 11.768310207409813, 'stdP': 192.79740339753218, 'stdRNA': 6.933396613857465, 'minP': 342.49949120254036, 'minRNA': 15.636285190573679, 'maxP': 796.5365881305912, 'maxRNA': 4.996731346683865}

summary_statistics = {
            k + species: lambda x: (s(x) - stat_mean[k + species])/stat_std[k + species]
            for k, s in summary_statistics_raw.items()
            for species in ['P', 'RNA']
        }

lst_square = lambda x,y: (
        (np.array([summary_statistics[k + species](x[species])
                   for k in summary_statistics_raw.keys()
                   for species in ['P', 'RNA']])
        - (np.array([summary_statistics[k + species](y[species])
                   for k in summary_statistics_raw.keys()
                   for species in ['P', 'RNA']]))
    )**2).sum()**0.5


with open(sys.argv[1], 'r') as f:
    data = json.load(f)
    y_observed = {
            'P':np.array(data[1]['P'])[:, BURNIN:].flatten().tolist(),
            'RNA':np.array(data[1]['RNA'])[:, BURNIN:].flatten().tolist(),
            }

abc = pyabc.ABCSMC(models, priors, lst_square)

db_path = "sqlite:///" + sys.argv[2].replace('(', '_').replace(')', '_').replace(',', '_')
history = abc.new(db_path, y_observed)
history = abc.run(minimum_epsilon=EPSILON, max_nr_populations=MAXPOP)
