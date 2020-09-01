from WMM import *

def test_generate_partitions():
    for param, result in [
            ((1, 1), [1]),
            ((10, 1), [10]),
            ((10, 2), [5, 5]),
            ((10, 3), [4, 3, 3]),
            ((3, 5), [1, 1, 1, 0, 0]),
            ]:
        assert generate_partitions(*param) == result

def test_run():
    """ Test no errors are raised

    """
    parameters = {
        'gamma': 1.,
        'mu': 1.,
        'kappa': 1.,
        'k_a': 1.,
        'k_d': 1.,
        }
    model = Cell(parameters)
    results = model.run({'Gf':1}, 100, 5, n_points=100)

def test_results_many_traj():
    """ Test results are collected properly

    """
    parameters = {
        'gamma': 0.,
        'mu': 0.,
        'kappa': 0.,
        'k_a': 0.,
        'k_d': 0.,
        }

    model = Cell(parameters)

    initial_state = {'Gf':1, 'Gb':2, 'RNA':3, 'P':4}

    n_traj = 5
    results = model.run(
            initial_state,
            time_stop=1, n_trajectories=n_traj, n_points=2)

    for species, n in initial_state.items():
        for i in range(n_traj):
            assert results[species][i][1] == n

def test_results_single_traj():
    """ Test results are collected properly

    """
    parameters = {
        'gamma': 0.,
        'mu': 0.,
        'kappa': 0.,
        'k_a': 0.,
        'k_d': 0.,
        }

    model = Cell(parameters)

    initial_state = {'Gf':1, 'Gb':2, 'RNA':3, 'P':4}

    results, final_state = model.run(
            initial_state,
            time_stop=1, n_trajectories=1, n_points=2)

    assert initial_state == final_state

    for species, n in initial_state.items():
        assert results[species][-1] == n

def deprecated_test_parallel():
    """ Test parallel run using dask

    """
    from dask.distributed import Client

    parameters = {
        'gamma': 1.,
        'mu': 1.,
        'kappa': 1.,
        'k_a': 1.,
        'k_d': 1.,
        'k_nc': 1.,
        'k_cn': 1.,
        }
    model = Cell(parameters)

    results = model.run({'Gf':1}, 100, 5, n_points=100, client=Client())

    assert len(list(results.values())[0]) == 5
