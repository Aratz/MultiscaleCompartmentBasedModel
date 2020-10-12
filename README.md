# A multiscale compartment-based model of stochastic gene regulatory networks using hitting-time analysis
Spatial stochastic models of single cell kinetics are capable of capturing both
fluctuations in molecular numbers and the spatial dependencies of the key steps
of intracellular regulatory networks. The spatial stochastic model is
straightforward to formulate and can be simulated using existing software, but
due to its computationally cost it quickly becomes prohibitively expensive.
This limits it use in applications that requires repeated simulation of the
model, such as when embedded in multicellular simulations, for parameter
inference, and in robustness analysis, model exploration and model checking. We
here propose a multiscale model where a compartment-based model approximates a
detailed spatial stochastic model. The compartment model is constructed via a
first-exit times analysis on the spatial model, thus capturing spatial aspects
of the fine-grained simulations. We apply the approach to a model of negative
feedback gene regulation and evaluate the approximation accuracy over a wide
range of parameters, assessing the situations in which a detailed spatial
representation can be replaced by the computationally much cheaper compartment
multiscale model.

# How to use this repository
This repository is organized as follows:
- The scripts used to generate the data are gathered in the `src` folder.
- All the data necessary to generate the figures is compressed in the `zip` files
from the `data` directory and must be extracted before the figures can be reproduced.
- All the figures from the manuscript can be regenerated using the Jupyter
notebooks.

All the computations were performed
on resources provided by SNIC through Uppsala Multidisciplinary Center for
Advanced Computational Science (UPPMAX) under Project SNIC 2019/8-227.
