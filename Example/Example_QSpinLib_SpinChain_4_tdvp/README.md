# QSpinLib Example\: Spin Chain ground state dynamics
We offer here an example of density matrix renormalization group (DMRG). \
Run **RunQSpinLib_DMRG.m** firstly to obtain the groud-state MPS. \
Run **RunQSpinLib_tdvp.m** secondly to calculate the ground-state spin dynamics.

## Hamiltonian ##
$H=\sum_{i=1}^{L-1} J_{xy} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y) + J_z S_i^z S_{i+1}^z - \sum_{i=1}^L (h_x S_i^x + h_y S_i^y + h_z S_i^z)$

## Results ##
### Results for ground state ###
* Ground-state energy $E_0$ ```Rslt.E0```
* Ground-state per-site magnetization $M$ ```Rslt.M```
* Ground-state MPS ```Rslt.T```
* Hamiltonian MPO ```Rslt.H```

### Results for real-time evolution ###
*


## Script and parameter settings ##
**Note: We have omitted the parameters that were mentioned earlier.**

### RunQSpinLib_DMRG.m ###

* Set the many-body solver as DMRG and start calculation.
```matlab
save(['MPS/', FileName], 'Para', 'Rslt')
```

### ImportMBSolverPara.m ###
* Save ground-state MPS and Hamiltonian MPO.
```matlab
Para.saveMPS = true;
```
