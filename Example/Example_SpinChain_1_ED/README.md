# Example\: Spin Chain ED
We offer here an example of exact diagonalization (ED). \
Please run **RunQSpinLib_ED.m** to start.

## Hamiltonian ##
$$H=\sum_{i=1}^{L-1} J_{xy} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y) + J_z S_i^z S_{i+1}^z - \sum_{i=1}^L (h_x S_i^x + h_y S_i^y + h_z S_i^z)$$

## Results ##
* Matrix of the Hamiltonians $H$ ```Rslt.H```
* Temperature $T$ ```Rslt.T```
* Inverse temperature $\beta \equiv \frac{1}{T}$ ```Rslt.beta```
* Per-site free energy $f = -\frac{1}{L \beta} log(\mathcal{Z})$ with $\mathcal{Z} \equiv {\rm Tr}[e^{-\beta H}]$ ```Rslt.F```
* Per-site energy $e = \frac{1}{L \mathcal{Z}} {\rm Tr}[e^{-\beta H} H]$ ```Rslt.En```
* Specific heat $C_{\rm m} = \frac{\partial E}{\partial T}$ ```Rslt.Cm```
* Magnetization $M = -\frac{\partial f}{\partial h}$ with $h = \sqrt{h_x^2 + h_y^2 + h_z^2}$ ```Rslt.M```
  
## Script and parameter settings ##

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

* Set the calculated thermodynamic quantity as specific heat.
```matlab
Para.ThDQ = 'Cm';
```
For the magnetic susceptibility $M/h$, please set
```matlab
Para.ThDQ = 'Chi';
```

* Set the many-body solver as exact diagonalization.
```matlab
Para.ManyBodySolver = 'ED';
```
For constructing the sparse matrix of Hamiltonian, please set
```matlab
Para.ManyBodySolver = 'ED_SM';
```
Note: the magnetic field is ignored in ```'ED_SM'```

* Set the output file to the standard output.
```matlab
Para.fileID = 1;
```

* Import the parameter of many-body solver. (see below)
```matlab
Para = ImportMBSolverPara(Para);
```

* Obtain the result and display the result.
```matlab
switch Para.ManyBodySolver
    case 'ED'
        Rslt = GetEDRslt(Para);
    case 'ED_SM'
        H = GetEDSMRslt(Para);
end

if strcmp(Para.ManyBodySolver, 'ED')
    switch Para.ThDQ
        case 'Cm'
            semilogx(Rslt.T, Rslt.Cm, '-o', 'linewidth', 2);
            xlabel('$T$', 'Interpreter', 'latex')
            ylabel('$C_{\rm m}$', 'Interpreter', 'latex')
            set(gca, 'XColor', 'k', 'YColor', 'k', 'fontsize', 20, 'fontname', 'times new roman', 'linewidth', 1.5)
        case 'Chi'
            semilogx(Rslt.T, Rslt.M/norm(Para.Field.h), '-o', 'linewidth', 2);
            xlabel('$T$', 'Interpreter', 'latex')
            ylabel('$\chi$', 'Interpreter', 'latex')
            set(gca, 'XColor', 'k', 'YColor', 'k', 'fontsize', 20, 'fontname', 'times new roman', 'linewidth', 1.5)
    end
end
```

### ImportMBSolverPara.m ###
Set the inverse temperature $\beta$ for ED.



