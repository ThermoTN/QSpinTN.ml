function UnitCon = GetUnitCon(para)
%GETUNITCON Package physical conversion, preserving native susceptibility.
u = QSpinThermalUnits();
ES = para.Model.ES;
UnitCon.T_con = ES;
UnitCon.Cm_con = u.R_J_per_mol_K;
UnitCon.h_con = ES*u.kB_J_per_K/u.muB_J_per_T;
% Native Chi is SI m^3/mol multiplied by 10^6, NOT CGS emu/mol/Oe.
UnitCon.Chi_con = u.NA_per_mol*u.mu0_SI*u.muB_J_per_T^2/ ...
    (ES*u.kB_J_per_K)*1e6;
end
