# QSpinLib Example\: Spin Chain ED
We offer here an example of exact diagonalization. \
Please run **RunQSpinLib_ED.m** to start.

## Hamiltonian ##
$H=\sum_{i=1}^{L-1} J_{xy} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y) + J_z S_i^z S_{i+1}^z - \sum_{i=1}^L (h_x S_i^x + h_y S_i^y + h_z S_i^z)$

## Workflow ##

### RunQSpinLib_ED.m ###
* Initialize the working environment.
```matlab
clear; clc;
maxNumCompThreads(4);
addpath(genpath('../../ManyBodySolver/'));
addpath(genpath('../../SpinModel/'));
addpath(genpath('../../svd_lapack_interface/'));
```
* Set model.
```matlab
Para.IntrcMap_Name = 'IntrcMap_XXZtest';
```
* Set coupling strength $J_{xy}=1$ and $J_z=1$.
```matlab
Para.Model.Jxy = 1; 
Para.Model.Jz = 1;
```

* Set the dimension of local Hilbert space. For spin-1/2 system, we have $d=2S+1=2$.
```matlab
Para.d = 2;
```
* Set lattice geometry with $L=12$ and open boundary condition.
```matlab
Para.Geo.L = 12;
Para.Geo.BC = 'OBC';
Para.L = Para.Geo.L;
```
For periodic boundary condition, plase set
```matlab
Para.Geo.BC = 'PBC';
```
Thus a term $J_{xy} (S_1^x S_{L}^x + S_1^y S_{L}^y) + J_z S_1^z S_{L}^z$ will be introduced.

* Set magnetic field as $h_x=0$, $h_y=0$, and $h_z=0.1$.
```matlab
Para.Field.h = [0,0,0.1];
```
For arbitrary magnetic field, please set
```matlab
Para.Field.h = [hx, hy, hz];
```

*
