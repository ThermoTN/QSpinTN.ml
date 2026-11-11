function physical = ConvertThermalObservables(raw, field_T, g)
%CONVERTTHERMALOBSERVABLES Natural per-spin observables to explicit units.
% raw.C_natural is C/kB per spin. raw.spin_projection is <S> along g.*B.
% Chi is finite-field M/H for a single nonzero principal-axis field.
% For zero/oblique fields Chi is NaN, never a differential susceptibility.
validateattributes(field_T,{'numeric'},{'vector','numel',3,'real','finite'});
validateattributes(g,{'numeric'},{'vector','numel',3,'real','finite','positive'});
u = QSpinThermalUnits();
physical = struct('Cm_J_per_mol_K',u.R_J_per_mol_K*raw.C_natural, ...
    'M_muB_per_spin',nan(size(raw.spin_projection)), ...
    'Chi_emu_per_mol_Oe',nan(size(raw.spin_projection)), ...
    'Chi_qspintn_native',nan(size(raw.spin_projection)), 'units',u);
axis = find(field_T~=0);
if numel(axis)==1
    physical.M_muB_per_spin = g(axis)*sign(field_T(axis))*raw.spin_projection;
    physical.Chi_emu_per_mol_Oe = u.NA_muB_emu_per_mol* ...
        physical.M_muB_per_spin/(u.Oe_per_T*field_T(axis));
    physical.Chi_qspintn_native = physical.Chi_emu_per_mol_Oe/u.native_chi_to_cgs;
elseif isempty(axis)
    physical.M_muB_per_spin = zeros(size(raw.spin_projection));
end
end
