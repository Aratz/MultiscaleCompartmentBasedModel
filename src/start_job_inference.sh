#!/bin/bash
sbatch launch_inference.sh data/smoldyn_data\(0.8000000000000007,0.8000000000000003,0.0\).json inference/inference\(0.8000000000000007,0.8000000000000003,0.0\).db
sbatch launch_inference.sh data/smoldyn_data\(-0.7999999999999998,0.8000000000000003,0.0\).json inference/inference\(-0.7999999999999998,0.8000000000000003,0.0\).db
sbatch launch_inference.sh data/smoldyn_data\(-4.8,3.2,0.0\).json inference\(-4.8,3.2,0.0\).db
