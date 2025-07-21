# Example\: Spin Chain XTRG
We offer here an example of exponential thermal tensor network (XTRG). \
**Note: Please read the article "Chen, et al. Phys. Rev. X. 8, 031082" before run the script.**\
Please run **RunQSpinLib_XTRG.m** to start. \

## Hamiltonian ##
$$H=\sum_{i=1}^{L-1} J_{xy} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y) + J_z S_i^z S_{i+1}^z - \sum_{i=1}^L (h_x S_i^x + h_y S_i^y + h_z S_i^z)$$

## Algorithm ##
We start from the density matrix at high temperature $\rho(\tau_0) = e^{-\tau_0 H}$ with $T=1/\tau_0$ and $\tau_0 << 1$.
Then we double the inverse temperature, i.e. $\rho(2\beta) = \rho(\beta)\rho^\dagger(\beta)$ and obtain the density matrix at 
$2 \tau_0,~4 \tau_0, ..., 2^n \tau_0$.


## Results ##
* Inverse temperature $\beta \equiv \frac{1}{T}$ ```Rslt.beta```
* Per-site $\frac{1}{L}log(\mathcal{Z})$ with $\mathcal{Z} \equiv {\rm Tr}[e^{-\beta H}]$ ```Rslt.LnZ```
* Per-site energy $e = \frac{1}{L \mathcal{Z}} {\rm Tr}[e^{-\beta H} H]$ ```Rslt.En```
* Magnetization $M = -\frac{\partial f}{\partial h}$ with $h = \sqrt{h_x^2 + h_y^2 + h_z^2}$ ```Rslt.M```
* Truncation error of tanTRG ```Rslt.TrunErr```
* Entanglement entropy of MPO ```Rslt.EE```
  
## Script and parameter settings ##
**Note: We have omitted the parameters that were mentioned earlier.**

### RunQSpinLib_tanTRG.m ###
* Set the many-body solver as XTRG and start calculation.
```matlab
Para.ThDQ = 'Cm';
Para.ManyBodySolver = 'XTRG';
Para.fileID = 1;
Para = ImportMBSolverPara(Para);

Rslt = GetXTRGRslt(Para);
```

### ImportMBSolvePara.m ###

* Set SETTN temperature.
```matlab
Para.tau0 = 0.00025;
```

* Set maximum iteration times of XTRG.
```matlab
Para.It = 20;
```
 
* Set bond dimension of SETTN.
```matlab
Para.MCrit = 100;
```

* Set bond dimension of XTRG at each iteration.
```matlab
Para.D_list = 1:1:Para.It;
Para.D_list(:) = 100;
Para.D_list(end:-1:end-4) = 100;
```

* Set maximum iteration times of varitional MPO product at each XTRG step.
```matlab
Para.sweep_step = 50;
```