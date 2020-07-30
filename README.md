## Introduction

This is a reference implementation of the following model:

  Pena, R.F.O., Zaks, M., Roque, A.C. (2018). 
  Spontaneous activity dynamics in random networks of spiking neurons with synaptic noise. 
  Journal of Computational Neuroscience, 45:1–28. doi 10.1007/s10827-018-0688-6

The very same code can has functions to reproduce results from the following papers  

  Tomov, P. , Pena, R.F. , Roque, A.C., Zaks, M.A. (2016).
  Mechanisms of self-sustained oscillatory states in hierarchical modular networks with mixtures of electrophysiological cell types. 
  Frontiers in Computational Neuroscience, 10:23. doi 10.3389/fncom.2016.00023
  
  Tomov, P. ,Pena, R.F. , Zaks, M.A. , Roque,A.C.(2014).
  Sustained oscillations, irregular firing, and chaotic dynamics in hierarchical modular networks with mixtures of electrophysiological cell types. 
  Frontiers in Computational Neuroscience, 8:103. doi 10.3389/fncom.2014.00103

## Platform information

**Platform:** Linux

**c++ (GCC):** 8.3.0

The network model is implemented using object-oriented programming in c++ .
For data processing and visualization we used standard functions available in Matlab.

## Code repository

This folder contains seven Python codes:
  *  **main.cpp:** Main script to run the simulations. Creates objects responsible for simulation and network structure.
  *  **simulation.cpp:** Class responsible to assemble network model structure and for running the simulation.
  *  **simulation.hpp:** Header file with definitions and macros.
  *  **izhi.cpp:** Class responsible for the neurons. Contain functions that take care of the integration. 
  *  **izhi.hpp:** Header file with definitions and macros.
  *  **alphasynapse.cpp:** Class responsible for synapses. It is used in simulation to create connections among neurons.
  *  **alphasynapse.hpp:** Header file with definitions and macros.

## Running the code

The main script used to simulate the network is from the terminal is:

  *  **run.sh:** Script where parameters of interest are changed. 

An example of how to run the scrip:

```
bash run.sh
```

Before running the scrip, the code must be compiled:

```
c++ *cpp -o code.out
```

After running the simulations, the code generates three files

raster.dat = matrix containing raster plot information. First column contains neuron index, then every line contains
 spike times. 

listsyn_M0_In_2_Ex2_3_Sim_121_NumSim_1.dat = Connectivity matrix

ex12in_M0_In_2_Ex2_3_Sim_121_NumSim_1.dat = Neuron type

We have included a Matlab code in order to create a figure of the raster plot. 
The default parameters should create a similar raster plot like the one in Fig.8 from the paper. 
It is expected that the figure won't be the same given its stochastic nature. 