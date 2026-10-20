# On the Practical Implementation of a Sequential Quadratic Programming Algorithm for Nonconvex Sum-of-squares Problems
[![Identifier](https://img.shields.io/badge/doi-10.18419%2Fdarus--5677-d45815.svg)](https://doi.org/10.18419/darus-5677)

This repository supplements the paper "On the Practical Implementation of a Sequential
Quadratic Programming Algorithm for Nonconvex
Sum-of-squares Problems" submitted to Mathematical Programming Computation. The preprint can be found on [arxiv](https://arxiv.org/pdf/2602.02394).
It contains the case benchmarks for different nonconvex sum-of-squares problems, with focus on analysis and control problems. The paper proposes a sequential quadratic programming like solver to tackle nonconvex sum-of-squares problems. 
This repository cotains the scripts to reproduce the benchmarks.

### General Requirements and Setup
The benchmarks itself does not not need to be installed. Only a [stable version](https://github.com/ifr-acso/casos/releases/latest) of CaΣoS, [MOSEK](https://www.mosek.com/downloads/) v11.1 is needed and CasADi v3.8.1 are required. 

  1. Download CasADi v3.8.1 and add it to your Matlab path.
  2. Download and install MOSEK, and add the solver to the Matlab path.
  3. Add the CaΣoS root folder (the one that contains the directory +casos) to your Matlab path.

If installed correctly, you can excute the corresponding files as outlined below for reproduction. The benchmarks in the  [paper](https://arxiv.org/pdf/2602.02394) were executed with MATLAB 2026a
on a computer with Windows 11 running on an AMD Ryzen 9 5950X 16-Core Processor with 3.40 GHz and 128 GB RAM.


### Folder Structure


```text
constraintViolationCheck_comparison.m   # Comparison from Section IV-A
00_ROA		                            # Folder containing the region-of-attraction cases
├── F18                                 # region of attraction for the 7D F18 aircraft
├── RobotArm                            # Contains the N-link robot arm example. 
01_Control_design                       # Folder containing for nonlinear systems
├── Aircraft_4D                         # Contains the control design for nonlinear longitudinal motion of an aircraft
├── Satellite_6D                        # Control desing for control-affine nonlinear dynamics of a satellite
02_Reachabilty                          # Folder contains two cases to compute an inner-approximation of the reachable set
├──Aircarft_4D                          # Inner-approximation of the reachable set for an aircraft lonitudinal motion
├──VanDerPol_2D                         # Inner-approximation of the reachable set for the Van-der-Pol Oscillator
03_CBF_CLF                              # Files for sequential SOS to compute a compatible pair of CBF-CLF
```
All folder contain the script for either the sequential SOS or coordinate-descent implementation. Except for the N-link robot arm we also provide screenshots showing the computational stats. Note that these might differ on other machines.

## Citation
Please cite the following paper, if you use the sequential algorithm.

> J. Olucak and T. Cunis , ‘On the Practical Implementation of a Sequential Quadratic Programming Algorithm for Nonconvex Sum-of-squares Problems’,  submitted to Mathematical Programming Computation, 2026. Pre-print:(https://arxiv.org/abs/2602.02394)


<details>

<summary>Bibtex entry</summary>

```bibtex
@misc{olucak2026practicalimplementationsequentialquadratic,
      title={On the Practical Implementation of a Sequential Quadratic Programming Algorithm for Nonconvex Sum-of-squares Problems}, 
      author={Jan Olucak and Torbjørn Cunis},
      year={2026},
      eprint={2602.02394},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2602.02394}, 
}
```



