function units = QSpinThermalUnits()
%QSPINTHERMALUNITS Shared physical constants and explicit output conventions.
% Constants preserve GetUnitCon's established numerical convention.
units = struct('schema_version','qspintn.thermal_observables.v1', ...
    'kB_J_per_K',0.08617*1.602e-22, 'muB_J_per_T',9.274e-24, ...
    'NA_per_mol',6.022e23, 'R_J_per_mol_K',8.31446, ...
    'mu0_SI',4*pi*1e-7, 'Oe_per_T',1e4, 'molar_basis','per mol spin');
units.muB_over_kB_K_per_T = units.muB_J_per_T/units.kB_J_per_K;
units.NA_muB_emu_per_mol = units.NA_per_mol*units.muB_J_per_T*1e3;
units.native_chi_to_cgs = 1/(4*pi);
units.native_chi_convention = 'qspintn_chi_si_cm3_v1';
end
