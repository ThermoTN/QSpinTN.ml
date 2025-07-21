# Example\: Triangular lattice ground-state dynamics
We offer here an example of ground-state spin dynamcis of triangular lattice quantum magnet.

## Hamiltonian ##
$$H=\sum_{\langle i,j \rangle} S_i^\alpha J_{i,j}^{\alpha \beta} S_j^\beta + \sum_{\langle \langle i,j \rangle \rangle } J_{2,xy} (S_i^x S_j^x + S_i^y S_j^y) + J_{2,z}S_i^z S_j^z - \sum_{i=1} (h_x S_i^x + h_y S_i^y + h_z S_i^z),$$
where $\langle i,j \rangle$ labels the nearest neighbor (NN) site, $\langle \langle i,j \rangle \rangle$ labels the next-nearest neighbor (NNN) site, and \
$$
J^{\alpha\beta} = 
\begin{pmatrix}
J_{xy}+2 J_{\rm PD} \cos\varphi
&-2J_{\rm PD} \sin\varphi
&-J_{\Gamma} \sin\varphi \\
-2J_{\rm PD} \sin\varphi
&J_{xy}-2J_{\rm PD} \cos\varphi
& J_{\Gamma} \cos\varphi\\
-J_{\Gamma} \sin\varphi
& J_{\Gamma} \cos\varphi
& J_z
\end{pmatrix}, \nonumber
$$

with $\varphi=\\{0,\frac{2\pi}{3},-\frac{2\pi}{3}\\}$ for three different types of NN bonds.

In **PlotSiteLocation.m**, we provide the lattice geometry, with the red, green, and blue dashed line corresponding to $\varphi=\\{0,\frac{2\pi}{3},-\frac{2\pi}{3}\\}$ respectively.


**Note: We have omitted the parameters that were mentioned earlier.**

Please see **Example_QSPinLib_TriARX_1_tanTRG** and **Example_QSpinLib_SpinChain_4_tdvp**.
