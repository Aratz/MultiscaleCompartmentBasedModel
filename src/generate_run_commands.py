"""
    Usage: `ls data/*_data*.json | xargs python3 generate_run_commands.py`
"""
import sys
import parse
import numpy as np


done = [
    parse.parse("data/{}_data({:f},{:f},{:f}).json", s).fixed
    for s in sys.argv[1:]] if len(sys.argv) > 1 else []
for solver in ['smoldyn', 'CBM', 'WMM']:
    params = set([
            (D, chi, 0)
            for D in np.linspace(-8, 4, num=16)
            for chi in np.linspace(-2, 4, num=16)
            if (solver, D, chi, 0) not in done]
        + [
            (D, 0, k_d)
            for D in np.linspace(-8, 4, num=16)
            for k_d in np.linspace(-3, 3, num=16)
            if (solver, D, 0, k_d) not in done])

    for D, chi, k_d in params:
        print(f"sbatch launch_{solver}.sh {D} {chi} {k_d}")
