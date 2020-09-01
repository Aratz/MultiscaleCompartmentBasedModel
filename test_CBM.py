import pytest
from CBM import *

@pytest.fixture
def parameters():
    return {
        'gamma': 1.,
        'mu': 1.,
        'kappa': 1.,
        'k_a': 1.,
        'k_d': 1.,
        'k_nc': 1.,
        'k_cn': 1.,
        }

@pytest.fixture
def null_parameters():
    return {
        'gamma': 0.,
        'mu': 0.,
        'kappa': 0.,
        'k_a': 0.,
        'k_d': 0.,
        'k_nc': 0.,
        'k_cn': 0.,
        }

def test_generate_partitions():
    for param, result in [
            ((1, 1), [1]),
            ((10, 1), [10]),
            ((10, 2), [5, 5]),
            ((10, 3), [4, 3, 3]),
            ((3, 5), [1, 1, 1, 0, 0]),
            ]:
        assert generate_partitions(*param) == result

def test_run(parameters):
    """ Test no errors are raised

    """
    model = Cell(parameters)
    results = model.run({'Gf':1}, 100, 5, n_points=100)

def test_results_many_traj(null_parameters):
    """ Test results are collected properly

    """
    model = Cell(null_parameters)

    initial_state = {
            'Gf':1, 'Gb':2, 'RNAnuc':3,
            'RNAcyt': 4, 'Pnuc':5, 'Pcyt':6}

    n_traj = 5
    results = model.run(initial_state,
            time_stop=1, n_trajectories=n_traj, n_points=2)

    for species, n in initial_state.items():
        for i in range(n_traj):
            assert results[species][i][1] == n

def test_results_single_traj(null_parameters):
    """ Test results are collected properly

    """
    model = Cell(null_parameters)

    initial_state = {
            'Gf':1, 'Gb':2, 'RNAnuc':3,
            'RNAcyt': 4, 'Pnuc':5, 'Pcyt':6}

    results, final_state = model.run(initial_state,
            time_stop=1, n_trajectories=1, n_points=2)

    assert initial_state == final_state

    for species, n in initial_state.items():
        assert results[species][-1] == n


def deprecated_test_parallel(parameters):
    """ Test parallel run using dask

    """
    import time
    from dask.distributed import LocalCluster, Client

    cluster = LocalCluster(processes=False)
    client = Client(cluster)

    model = Cell(parameters)

    s_start = time.time()
    s_results = model.run({'Gf':1}, 1000, 1024, n_points=100)
    s_end = time.time()

    p_start = time.time()
    p_results = model.run({'Gf':1}, 1000, 1024, n_points=100, client=client)
    p_end = time.time()

    assert s_end - s_start > 0.9 * 4 * (p_end - p_start)

def test_npoints(parameters):
    """ Test number of samples per trajectory

    """
    import numpy as np

    n_points = 100

    # Test with one trajectory
    n_traj = 1
    for n_traj in [1, 5, 10]:
        model = Cell(parameters)
        results = model.run({'Gf':1}, 100, n_traj, n_points=n_points)
        history = results[0] if n_traj == 1 else results
        for species in history:
            assert np.array(history[species]).shape == (
                (n_points + 1,) if n_traj == 1 else (n_traj, n_points + 1))
