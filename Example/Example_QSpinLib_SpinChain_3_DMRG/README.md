# QSpinLib Example\: Spin Chain DMRG
We offer here an example of density matrix renormalization group (DMRG). \
Please run **RunQSpinLib_DMRG.m** to start. \

## Hamiltonian ##
$H=\sum_{i=1}^{L-1} J_{xy} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y) + J_z S_i^z S_{i+1}^z - \sum_{i=1}^L (h_x S_i^x + h_y S_i^y + h_z S_i^z)$

## Results ##
* Ground-state energy $E_0$ ```Rslt.E0```
* Ground-state per-site magnetization $M$ ```Rslt.M```
  
## Script and parameter settings ##
**Note: We have omitted the parameters that were mentioned earlier.**

### RunQSpinLib_DMRG.m ###

* Set the many-body solver as DMRG and start calculation.
```matlab
Para.ManyBodySolver = 'DMRG';
Para.fileID = 1;
Para = ImportMBSolverPara(Para);

Rslt = GettanTRGRslt(Para);
```
Note: Please do not modify ```Para.ThDQ```.

### ImportMBSolverPara.m ###
* Set the bond dimesnion of MPS.
```matlab
Para.MCrit = 100;
```

* Set the dimension of Lanczos space.
```matlab
Para.DK = 20;
```

* Set the tolerance of Lanczos process.
```matlab
Para.tol = 1e-12;
```

* Set the maximum number of DMRG iterations.
```matlab
Para.DMRGStepMax = 100;
```
