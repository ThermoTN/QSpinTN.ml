# QSpinLib Example\: Spin Chain tanTRG with DM interaction
We offer here an example of tangent space tensor renormalization group (tanTRG) for 1D spin chain with DM interactions.
The main purpose of this example is to demonstrate how to construct the Hamiltonian.

## Hamiltonian ##
$$H=\sum_{i=1}^{L-1} J_{xy} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y) + J_z S_i^z S_{i+1}^z + D\cdot(S_i\times S_{i+1})- \sum_{i=1}^L (h_x S_i^x + h_y S_i^y + h_z S_i^z)$$

## Script and parameter settings ##
**Note: We have omitted the parameters that were mentioned earlier.**

### RunQSpinLib_tanTRG.m ###
* Set the interaction map and the coupling strenght.
```matlab
Para.IntrcMap_Name = 'IntrcMap_1DDMChain';
Para.Model.Jxy = 1;
Para.Model.Jz = 1;
Para.Model.DM = [1,1,1];
```

### IntrcMap_1DDMChain.m ###
* Construct the interaction map of the Hamiltonian.
```matlab
[ Intr ] = IntrcMap_1DDMChain( Para )
```
```Intr``` is a **struch** with $n$ terms, with $n$ the number of interactions. \
For each term (e.g. the $i$th term), it has three fields, i.e.
```matlab
Intr(i).site = [l1, l2, ..., lm];
Intr(i).operator = {Op1, Op2, ..., Opm};
Intr(i).J = J;
```
where $1 \leq l_1 < l_2 <... < l_m \leq L$ labels the site of the corresponding operators ($d \times d$ matrix) ${\rm Op}_1, {\rm Op}_2, ..., {\rm Op}_m$ 
with the coupling strenght $J$.
For example:
* $0.5 S_1^x S_2^z$
```matlab
Intr(i).site = [1, 2];
Intr(i).operator = {[0,0.5;0.5,0], [-0.5,0;0,0.5]};
Intr(i).J = 0.5;
```
* $-0.2 S_5^x$
```matlab
Intr(i).site = [5];
Intr(i).operator = {[0,0.5;0.5,0]};
Intr(i).J = -0.2;
```
* $-(S^z)^2$ (spin-1 system)
```matlab
Intr(i).site = [5];
Intr(i).operator = {[1,0,0;0,0,0;0,0,1]};
Intr(i).J = -1;
```
* $S_1^x S_2^x S_5^x S_6^x$
```matlab
Intr(i).site = [1, 2, 5, 6];
Intr(i).operator = {[0,0.5;0.5,0], [0,0.5;0.5,0], [0,0.5;0.5,0], [0,0.5;0.5,0]};
Intr(i).J = 1;
```
