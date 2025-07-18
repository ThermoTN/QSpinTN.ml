# QSpinLib Example\: Spin Chain tanTRG
We offer here an example of tangent space tensor renormalization group (tanTRG). \
Please run **RunQSpinLib_tanTRG.m** to start. \
We have already provided the result in
```matlab
Rslt/Model=XXZtest-Jxy=1-Jz=1-h=[0 0 0.1]-L=12-BC=OBC-D=200.mat
```
To save time, you can directly run **Benchmark.m**.


## Hamiltonian ##
$H=\sum_{i=1}^{L-1} J_{xy} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y) + J_z S_i^z S_{i+1}^z - \sum_{i=1}^L (h_x S_i^x + h_y S_i^y + h_z S_i^z)$

## Results ##
* Inverse temperature $\beta \equiv \frac{1}{T}$ ```Rslt.beta```
* Per-site $log(\mathcal{Z})$ with $\mathcal{Z} \equiv {\rm Tr}[e^{-\beta H}]$ ```Rslt.LnZ```
* Per-site energy $E = \frac{1}{L \mathcal{Z}} {\rm Tr}[e^{-\beta H} H]$ ```Rslt.En```
* Magnetization $M = -\frac{\partial f}{\partial h}$ with $h = \sqrt{h_x^2 + h_y^2 + h_z^2}$ ```Rslt.M```
* Truncation error of tanTRG ```Rslt.TrunErr```
* Entanglement entropy of MPO ```Rslt.EE```
  
## Script and parameter settings ##

### RunQSpinLib_tanTRG.m ###
**Note: We have omitted the parameters that were mentioned earlier.**
* Set the many-body solver as tanTRG and start calculation.
```matlab
Para.ThDQ = 'Cm';
Para.ManyBodySolver = 'tanTRG';
Para.fileID = 1;
Para = ImportMBSolverPara(Para);

Rslt = GettanTRGRslt(Para);
```
Note: Please do not modify ```Para.ThDQ```.

### ImportMBSolverPara.m ###
* Set the temperature $\tau_0$ for SETTN initialization.
```matlab
Para.tau0 = 2^-15;
```

* Set the temperature step of tanTRG.
```matlab
Para.beta_max = 32;
[Para.tau_step, Para.beta_c ] = GetStepFunc(Para.tau0, [1, 0.25], Para.beta_max);
```
For the present case, the temperature will exponentially (double) increase to 1, and then linearly increase to 32 (with step $4\times0.25=1$). 
Please see ```Rslt.beta``` for more details. 

**For advanced users:**
We provide a method for arbitrarily designing the step size as
```matlab
[Para.tau_step, Para.beta_c ] = GetStepFunc(Para.tau0, [B, step1, step2, ..., stepN], [beta1, beta2, ..., beta_max]);
```
Thus the temperature will double increase to **B**, and then linearly increase to **beta1** with step **4*step1**, 
and then linearly increase to **beta2** with step **4*step2** and so on.
