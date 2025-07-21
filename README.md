# QSpinTN.ml (v0.0)
## Introduction ##
QSpinTN.ml is a package for the many-body simulations of the spin system (both ground state and finite-temperature) and an advanced tool for thermal data analysis of magnetic quantum materials.

## Quantum Many-Body Sovlers ##
* Exact diagonalization (ED, as a high-*T* solver);
* Density matrix renormalization group (DMRG, as a ground state solver for quasi-1D and 2D magnets);
* Linearized tensor renormalization group (LTRG, as a low-*T* solver for 1D spin chain materials, etc);
* Exponential tensor renormalization group (XTRG, as a low-*T* solver for quasi-1D and 2D magnets, etc).
* Tangent space tensor renormalization group (tanTRG, as a low-*T* solver for quasi-1D and 2D magnets, etc).

## Try Your First QSpinTN.ml Program ##
Here we provide some examples in **Example/**. Feel free to give it a try!

## Maintainer ##
* Yuan Gao, Beihang University\
  e-mail: yuangao@buaa.edu.cn
* Wei Li, ITP-CAS\
  e-mail: w.li@itp.ac.cn
  
## Acknowledgement ##
* Sizhuo Yu, CentraleSupélec
* Bin-Bin Chen, Beihang University
* Qiaoyi Li, ITP-CAS
* Dai-Wei Qu, UCAS
* Han Li, Beijing Jiaotong University
* Enze Lv, ITP-CAS
## Reference ##

```bib
@article{QMagen2021,
   title = {Learning the Effective Spin {H}amiltonian of a Quantum Magnet},
   author = {Sizhuo Yu and Yuan Gao and Bin-Bin Chen and Wei Li},
   publisher = {Chin. Phys. Lett.},
   year = {2021},
   journal = {Chin. Phys. Lett.},
   volume = {38}, 
   number = {9},
   eid = {097502},
   pages = {097502},
   url = {http://cpl.iphy.ac.cn/EN/abstract/article_115992.shtml},
   doi = {10.1088/0256-307X/38/9/097502}
}   
```
```bib
@article{LTRG2011,
  title = {Linearized Tensor Renormalization Group Algorithm for the Calculation of Thermodynamic Properties of Quantum Lattice Models},
  author = {Li, W. and Ran, S.-J. and Gong, S.-S. and Zhao, Y. and Xi, B. and Ye, F. and Su, G.},
  journal = {Phys. Rev. Lett.},
  volume = {106},
  issue = {12},
  pages = {127202},
  numpages = {4},
  year = {2011},
  month = {Mar},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.106.127202},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.106.127202}
}
```
```bib
@article{BilayerLTRG2017,
	Author = {Dong, Y.-L. and Chen, L. and Liu, Y.-J. and Li, W.},
	Doi = {10.1103/PhysRevB.95.144428},
	Issue = {14},
	Journal = {Phys. Rev. B},
	Month = {Apr},
	Numpages = {10},
	Pages = {144428},
	Publisher = {American Physical Society},
	Title = {Bilayer linearized tensor renormalization group approach for thermal tensor networks},
	Url = {https://link.aps.org/doi/10.1103/PhysRevB.95.144428},
	Volume = {95},
	Year = {2017},
	Bdsk-Url-1 = {https://link.aps.org/doi/10.1103/PhysRevB.95.144428},
	Bdsk-Url-2 = {http://dx.doi.org/10.1103/PhysRevB.95.144428}
}
```
```bib
@article{SETTN2017,
	Author = {Chen, B.-B. and Liu, Y.-J. and Chen, Z. and Li, W.},
	Doi = {10.1103/PhysRevB.95.161104},
	Issue = {16},
	Journal = {Phys. Rev. B},
	Month = {Apr},
	Numpages = {5},
	Pages = {161104(R)},
	Publisher = {American Physical Society},
	Title = {Series-expansion thermal tensor network approach for quantum lattice models},
	Url = {https://link.aps.org/doi/10.1103/PhysRevB.95.161104},
	Volume = {95},
	Year = {2017},
	Bdsk-Url-1 = {https://link.aps.org/doi/10.1103/PhysRevB.95.161104},
	Bdsk-Url-2 = {http://dx.doi.org/10.1103/PhysRevB.95.161104}
}
```
```bib
@Article{XTRG2018,
  title = {Exponential Thermal Tensor Network Approach for Quantum Lattice Models},
  author = {Chen, B.-B. and Chen, L. and Chen, Z. and Li, W. and Weichselbaum, A.},
  journal = {Phys. Rev. X},
  volume = {8},
  issue = {3},
  pages = {031082},
  numpages = {29},
  year = {2018},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevX.8.031082},
  url = {https://link.aps.org/doi/10.1103/PhysRevX.8.031082}
}
```
```bib
@article{XTRG2019,
  title = {Thermal tensor renormalization group simulations of square-lattice quantum spin models},
  author = {Li, H. and Chen, B.-B. and Chen, Z. and von Delft, J. and Weichselbaum, A. and Li, W.},
  journal = {Phys. Rev. B},
  volume = {100},
  issue = {4},
  pages = {045110},
  numpages = {17},
  year = {2019},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevB.100.045110},
  url = {https://link.aps.org/doi/10.1103/PhysRevB.100.045110}
}
```
```bib
@article{tanTRG2023,
  title = {Tangent Space Approach for Thermal Tensor Network Simulations of the 2D {H}ubbard Model},
  author = {Li, Qiaoyi and Gao, Yuan and He, Yuan-Yao and Qi, Yang and Chen, Bin-Bin and Li, Wei},
  journal = {Phys. Rev. Lett.},
  volume = {130},
  issue = {22},
  pages = {226502},
  numpages = {8},
  year = {2023},
  month = {Jun},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.130.226502},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.130.226502}
}
```
```bib
@article{TDVP2011,
  title = {Time-Dependent Variational Principle for Quantum Lattices},
  author = {Haegeman, Jutho and Cirac, J. Ignacio and Osborne, Tobias J. and Pi\ifmmode \check{z}\else \v{z}\fi{}orn, Iztok and Verschelde, Henri and Verstraete, Frank},
  journal = {Phys. Rev. Lett.},
  volume = {107},
  issue = {7},
  pages = {070601},
  numpages = {5},
  year = {2011},
  month = {Aug},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.107.070601},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.107.070601}
}
```
```bib
@article{TDVP2016,
  title = {Unifying Time Evolution and Optimization with Matrix Product States},
  author = {Haegeman, Jutho and Lubich, Christian and Oseledets, Ivan and Vandereycken, Bart and Verstraete, Frank},
  journal = {Phys. Rev. B},
  volume = {94},
  issue = {16},
  pages = {165116},
  numpages = {10},
  year = {2016},
  month = {Oct},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevB.94.165116},
  url = {https://link.aps.org/doi/10.1103/PhysRevB.94.165116}
}
```
```bib
@article{DMRG1992,
  title = {Density matrix formulation for quantum renormalization groups},
  author = {White, Steven R.},
  journal = {Phys. Rev. Lett.},
  volume = {69},
  issue = {19},
  pages = {2863--2866},
  numpages = {0},
  year = {1992},
  month = {Nov},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.69.2863},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.69.2863}
}
```
