# Example\: 1D XY model at finite temperature
We offer here an example of tangent space tensor renormalization group (tanTRG). \
**Note: Please read the article "Li, et al. Phys. Rev. Lett. 130, 226502" before run the script.**\
Please run **RunQSpinLib_tanTRG.m** to start. \
We have already provided the result in
```matlab
% tanTRG result
Rslt/Model=XXZtest-Jxy=1-Jz=0-h=[0 0 0]-L=50-BC=OBC-D=50.mat
Rslt/Model=XXZtest-Jxy=1-Jz=0-h=[0 0 0]-L=50-BC=OBC-D=100.mat
Rslt/Model=XXZtest-Jxy=1-Jz=0-h=[0 0 0]-L=50-BC=OBC-D=200.mat
% XTRG result
Rslt/Model=XXZtest-Jxy=1-Jz=0-h=[0 0 0]-L=50-BC=OBC-D=50-XTRG.mat
Rslt/Model=XXZtest-Jxy=1-Jz=0-h=[0 0 0]-L=50-BC=OBC-D=100-XTRG.mat
```
To save time, you can directly run **Check_energy.m** and **Check_free_energy**.


## Hamiltonian ##
$$H=\sum_{i=1}^{L-1} J_{xy} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y)$$

## Results ##
* Inverse temperature $\beta \equiv \frac{1}{T}$ ```Rslt.beta```
* Per-site $\frac{1}{L}log(\mathcal{Z})$ with $\mathcal{Z} \equiv {\rm Tr}[e^{-\beta H}]$ ```Rslt.LnZ```
* Per-site energy $e = \frac{1}{L \mathcal{Z}} {\rm Tr}[e^{-\beta H} H]$ ```Rslt.En```
* Magnetization $M = -\frac{\partial f}{\partial h}$ with $h = \sqrt{h_x^2 + h_y^2 + h_z^2}$ ```Rslt.M```
* Truncation error of tanTRG ```Rslt.TrunErr```
* Entanglement entropy of MPO ```Rslt.EE```