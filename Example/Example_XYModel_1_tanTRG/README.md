# Example\: 1D XY model tanTRG
We offer here an example of tangent space tensor renormalization group (tanTRG). \
**Note: Please read the article "Li, et al. Phys. Rev. Lett. 130, 226502" before run the script.**\
Please run **RunQSpinLib_tanTRG.m** to start. \
We have already provided the result in
```matlab
Rslt/Model=XXZtest-Jxy=1-Jz=0-h=[0 0 0]-L=50-BC=OBC-D=50.mat
Rslt/Model=XXZtest-Jxy=1-Jz=0-h=[0 0 0]-L=50-BC=OBC-D=100.mat
```
To save time, you can directly run **Check_energy.m** and **Check_free_energy**.


## Hamiltonian ##
$$H=\sum_{i=1}^{L-1} J_{xy} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y)$$