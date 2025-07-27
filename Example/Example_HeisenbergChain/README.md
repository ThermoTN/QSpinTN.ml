# Example\: Ground-state dynamical spin structure factor of 1D Heisenberg Chain
We offer here an example of dynamical spin structure factor (DSS). \
**Note: Please read** ```Example_SpinChain_4_tdvp``` **firstly**.\
Run **RunQSpinLib_DMRG.m** firstly to obtain the groud-state MPS. \
Run **RunQSpinLib_tdvp.m** secondly to calculate the ground-state spin dynamics.
Run **PlotSkomega.m** to plot the DSS $S(k,\omega)$ of 1D Heisenberg Chain with lattice size $L=64$ and bond dimension $D=200$.

## Hamiltonian ##
$$H=\sum_{i=1}^{L-1} S_i^x S_{i+1}^x + S_i^y S_{i+1}^y + S_i^z S_{i+1}^z$$

## Dynamical Spin Structure factor ##

For the isotropic Heisenberg model, we have
$$S(k,\omega) = 3\int_{-\infty}^{\infty}dt \, e^{i\omega t} 
g(k,t), $$
with the correlation function $g(k,t)=\langle e^{iHt}S^z_{-k}e^{-iHt}S_k^z\rangle$ and $S_k^z = \frac{1}{\sqrt{L}}\sum_m^Le^{imk}S_m^z$. 

Note that $g^*(k,t) =\langle e^{-iHt}S^z_{-k}e^{iHt}S_k^z\rangle =g(k,-t)$, we have 
$$S(k,\omega)=3[\int_{0}^{\infty} dt\, e^{i\omega t}g(k,t) + \int_{0}^{\infty} dt\, e^{-i\omega t}g(k,-t)]=6\int_0^{\infty}dt\, {\rm Re}[e^{i\omega t}g(k,t)].$$

### Window function ### 
In practice, we have a maximum evolution time $t_{\rm max}$ in tdvp. Here we introduce a window function $W(t/t_{\rm max})$ which vanishes with $|t/t_{\rm max}|>1$.

Now we have 
$$S(k,\omega)\simeq 6\int_0^{\rm t_{\rm max}} dt\, {\rm Re}[e^{i \omega t} g(k,t)] W(t/t_{\rm max}).$$

Here we provide two type window function, the Parzen window function and the staircase function, demonstrated in **WindowFunction.m**. Notably, the staircase function may lead to sub peaks, which may confuse the researcher. 

### Energy resolution ###
The finite $t_{\rm max}$ leads to a finite energy resolution. To estimate the energy resolution, we consider the full width at half maximum of the Fourier transform of the window function. In practice, the energy resolution of Parzen window function is $\varepsilon\simeq 8/t_{\rm max}$. Please run **Energyresolution.m** to see the detail.
