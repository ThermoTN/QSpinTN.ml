# QSpinLib Example\: Spin Chain ground state dynamics
We offer here an example of density matrix renormalization group (DMRG). \
Run **RunQSpinLib_DMRG.m** firstly to obtain the groud-state MPS. \
Run **RunQSpinLib_tdvp.m** secondly to calculate the ground-state spin dynamics.

## Hamiltonian ##
$H=\sum_{i=1}^{L-1} J_{xy} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y) + J_z S_i^z S_{i+1}^z - \sum_{i=1}^L (h_x S_i^x + h_y S_i^y + h_z S_i^z)$

### Algorithm ###
We firstly obtain the matrix product state (MPS) of the ground state. And then consider the real-time correlation function \
$\langle e^{iHt}A^\dagger e^{-iHt} A \rangle$, \
where $ e^{-iHt} A|\psi \rangle$ is obtained by using the time-dependent varitional principle (tdvp).

## Results ##
### Results for ground state ###
* Ground-state energy $E_0$ ```Rslt.E0```
* Ground-state per-site magnetization $M$ ```Rslt.M```
* Ground-state MPS ```Rslt.T```
* Hamiltonian MPO ```Rslt.H```

### Results for real-time evolution ###
* Real time $t$ ```Rslt.time_list```
* Correlation function ```Rslt.Obs```
* Entanglement entropy ```Rslt.EE```

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

### RunQSpinLib_tdvp.m ###
* Load ground-state MPS and Hamiltonian MPO.
```matlab
load MPS/'Model=XXZtest-Jxy=1-Jz=1-h=[0 0 0.1]-L=10-BC=OBC-D=100.mat'
```
* Set quench operator $A$.
```matlab
Para.DY.k = pi;
Para.DY.QulocOp = [-0.5, 0; 0, 0.5];
```
For the present case, we take $A=\frac{1}{\sqrt{L}} \sum_{n=1}^{L} e^{ikn} S_n^z$ with $k=\pi$.

* Set tdvp parameters.

Time step $\tau=0.2$ (2$\times$Para.tau) 
``` matlab 
Para.tau = 0.1;
```
Maximum evolutionary time
```matlab
Para.tmax = 10;
```
* Obtain the interaction map of $A$.
```matlab
Para = InitKSpaceIntr(Para);
```

* Perform real-time evolution.
```matlab
GettdvpRslt(Para, Rslt.T, Rslt.H )
```

### InitKSpaceIntr.m ###
Generate the interaction map of the quench operator $A$.
```matlab
for it = 1:Para.L
    Intr(it).site = [it];
    Intr(it).operator = {Para.DY.QulocOp};
    Intr(it).J = exp(1i * Para.DY.k * it)/sqrt(Para.L);
end
```
For each iteration, we add an operator ``` exp(1i * Para.DY.k * it)/sqrt(Para.L) * Para.DY.QulocOp``` on site ```it```. 
