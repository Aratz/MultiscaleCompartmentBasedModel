import pytest
from smoldyn import *
from smoldyn import Mol, Cell


class Test_Mol:
    def test_smoldyn_input(self):
        assert Mol("RNA", 1.0, 2.0, 3.0).to_smoldyn_input() == "mol 1 RNA 1.0 2.0 3.0"
        assert Mol("P", 5, -10.0, 0.0).to_smoldyn_input() == "mol 1 P 5.0 -10.0 0.0"
        assert Mol("Gf", -0.2, -250, -1).to_smoldyn_input() == "mol 1 Gf -0.2 -250.0 -1.0"

    def test_repr(self):
        assert str(Mol("RNA", 1.0, 2.0, 3.0)) == "{'name': 'RNA', 'x': 1.0, 'y': 2.0, 'z': 3.0}"

def test_parse_last_state():
    raw_last_state = [
            "RNA(solution) -1 -4 6 191",
            "RNA(solution) 1 2 3 193",
            "Gf(solution) 0 0 0 192",
            ]
    last_state = parse_last_state(raw_last_state)
    assert Mol("RNA", -1.0, -4.0, 6.0) in last_state
    assert Mol("RNA", 1.0, 2.0, 3.0) in last_state
    assert Mol("Gf", 0.0, 0.0, 0.0) in last_state

def test_kolmogorov_distance():
    data1 = {
            'P': [1., 1., 1.],
            'RNA': [2., 2., 2.]}
    data2 = {
            'P': [3., 3., 3.],
            'RNA': [2., 2., 2.]}
    data3 = {
            'P': [3., 3., 3.],
            'RNA': [4., 4., 4.]}

    assert kolmogorov_distance(data1, data1) == 0.
    assert kolmogorov_distance(data1, data2) == 1.
    assert kolmogorov_distance(data1, data3) == 2.

class Test_Cell:
    @pytest.fixture
    def template(self):
        return """random_seed {seed} \n define K_A {k_a} \n define K_D {k_d} \n define MU {mu} \n define KAPPA {kappa} \n define GAMMA {gamma} \n define OUTPUT stdout \n \n dim 3 \n boundaries x -8 8 \n boundaries y -8 8 \n boundaries z -8 8 \n time_start 0 \n time_stop {time_stop} \n time_step {time_step} \n \n species Gf Gb RNA P \n difc RNA {diffusion} \n difc P {diffusion} \n display_size all(all) 0.02 \n color P(all) lightblue \n color RNA(all) navy \n color Gf(all) scarlet \n color Gb(all) darkred \n \n graphics opengl_good \n frame_thickness 0 \n \n start_surface membrane \n action both all reflect \n color both blue \n thickness 0.01 \n panel sphere 0 0 0 {cell_radius} 30 30 \n end_surface \n \n start_surface nucleus \n action front RNA reflect \n action back P reflect \n color both red \n thickness 0.01 \n panel sphere 0 0 0 {nucleus_radius} 30 30 \n end_surface \n \n start_compartment inside \n surface nucleus \n point 0 0 0\n end_compartment \n \n start_compartment cytoplasm \n surface membrane \n point 0 0 0\n compartment andnot inside \n end_compartment \n \n reaction transcription Gf -> Gf + RNA MU \n reaction compartment=cytoplasm translation RNA -> RNA + P KAPPA \n reaction binding Gf + P <-> Gb K_A K_D \n reaction compartment=cytoplasm degradation RNA|P -> 0 GAMMA \n reaction compartment=inside degradation2 RNA|P -> 0 GAMMA \n {initial_state} \n \n output_files OUTPUT \n cmd B molcountheader OUTPUT \n cmd N {freq_point} molcount OUTPUT \n cmd A echo OUTPUT "--Simulation ends--\\n" \n cmd A listmols OUTPUT \n \n end_file \n"""

    @pytest.fixture
    def parameters(self):
        return {
            'k_a': 0.001,
            'k_d': 0.1,
            'mu': 1,
            'kappa': 1,
            'gamma': 0.1,
            'diffusion': 0.002,
            'time_step': 0.01,
            'cell_radius': 1,
            'nucleus_radius': 0.3,
            'model': 'smoldyn',
            }


    def test_smoldyn(self, template, parameters):
        """ Test smoldyn is running without raising exceptions

        """

        initial_state = [Mol("Gf", 0.0, 0.0, 0.0)]

        model = Cell(template, parameters)
        model.run(initial_state, 10., n_trajectories=1)

    def test_smoldyn_multiple(self, template, parameters):
        """ Test smoldyn is running several trajectories without raising
        exceptions

        """

        initial_state = [Mol("Gf", 0.0, 0.0, 0.0)]

        model = Cell(template, parameters)
        model.run(initial_state, 10., n_trajectories=10)

    def test_smoldyn_run_trajectory(self, template, parameters):
        """ Test smoldyn runs and terminates in threshold mode

        """

        n_points = 500

        initial_state = [Mol("Gf", 0.0, 0.0, 0.0)]

        model = Cell(template, parameters)
        history, _ = run_trajectory(
                model.model,
                time_stop=10.,
                time_step=0.01,
                initial_state=initial_state,
                seed=0,
                n_points=n_points)

        assert len(history['P']) == n_points + 1
        #if it works, rerun with smaller tstep and retest

        history, _ = run_trajectory(
                model.model,
                time_stop=10.,
                time_step=0.01/4,
                initial_state=initial_state,
                seed=0,
                n_points=n_points)
        assert len(history['P']) == n_points + 1

    def test_smoldyn_threshold(self, template, parameters):
        """ Test smoldyn runs and terminates in threshold mode

        """

        n_points = 500

        initial_state = [Mol("Gf", 0.0, 0.0, 0.0)]

        model = Cell(template, parameters)
        data, time = model.run_threshold(initial_state, 10.,
                n_trajectories=10, threshold=0.25, n_points=n_points)

        for tstep in data:
            assert len(tstep['P'][0]) == n_points + 1

    def test_smoldyn_docker(self, template, parameters):
        """ Test smoldyn is running in docker without raising exceptions

        """

        initial_state = [Mol("Gf", 0.0, 0.0, 0.0)]

        model = Cell(template, parameters, docker='smoldyn')
        model.run(initial_state, 10., n_trajectories=1)

    def test_smoldyn_parallel(self, template, parameters):
        """ Test smoldyn parallel run using dask

        """
        import time
        from dask.distributed import LocalCluster, Client

        cluster = LocalCluster(processes=False)
        client = Client(cluster)

        initial_state = [Mol("Gf", 0.0, 0.0, 0.0)]

        model = Cell(template, parameters)

        s_start = time.time()
        model.run(initial_state, 100., n_trajectories=32)
        s_end = time.time()

        p_start = time.time()
        model.run(initial_state, 100., n_trajectories=32, client=client)
        p_end = time.time()

        assert s_end - s_start > 0.8 * 2 * (p_end - p_start)
        # Minimal speed up is a trade off here between garanteed performance
        # and time necessary to run the test

        client.close()
        cluster.close()

    def test_noseed(self, template, parameters):
        """ Test results are different when no seed is provided

        """
        initial_state = [Mol("Gf", 0.0, 0.0, 0.0)]

        model = Cell(template, parameters)

        # Run models one after the other to ensure timebased seed is different
        results = [model.run(initial_state, 10., n_trajectories=1)]
        results.append(model.run(initial_state, 10., n_trajectories=1))

        assert results[0] != results[1]


    def test_seed(self, template, parameters):
        """ Test results are the same when the same seed is provided

        """
        initial_state = [Mol("Gf", 0.0, 0.0, 0.0)]

        model = Cell(template, parameters)
        results = [
            model.run(initial_state, 10., n_trajectories=1, seed=0)
            for _ in range(2)]

        assert results[0] == results[1]

    def test_npoints(self, template, parameters):
        """ Test timeseries have the specified length

        """
        import numpy as np

        initial_state = [Mol("Gf", 0.0, 0.0, 0.0)]

        # Vary tstop
        n_points = 100

        model = Cell(template, parameters)
        for tstop in [1., 10., 100.]:
            results = model.run(
                    initial_state, tstop,
                    n_trajectories=1, n_points=n_points, seed=0)
            for species in results[0]:
                assert len(results[0][species]) == n_points + 1

        # Vary tstep
        n_points = 500
        tstop = 1000

        for tstep in [0.01, 0.02, 0.04]:
            parameters['time_step'] = tstep
            model = Cell(template, parameters)
            results = model.run(
                    initial_state, tstop,
                    n_trajectories=1, n_points=n_points, seed=0)
            for species in results[0]:
                assert len(results[0][species]) == n_points + 1

        # Test with many trajectories
        n_points = 500
        tstop = 1000
        n_traj = 10

        for tstep in [0.01, 0.02, 0.04]:
            parameters['time_step'] = tstep
            model = Cell(template, parameters)
            results = model.run(
                    initial_state, tstop,
                    n_trajectories=n_traj, n_points=n_points, seed=0)
            for species in results:
                assert np.array(results[species]).shape == (n_traj, n_points + 1)
