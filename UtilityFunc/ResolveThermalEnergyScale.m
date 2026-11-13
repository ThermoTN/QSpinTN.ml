function scale = ResolveThermalEnergyScale(coefficients_K, fields_T, g, g_axes)
%RESOLVETHERMALENERGYSCALE Common reference for one physical parameter point.
% coefficients_K contains actual Hamiltonian coefficients after expansion of
% dimensionless ratios and model factors, not MPO singular values or norms.
% Rows of fields_T are ALL fields evaluated at this point. Rows of g_axes
% are orthonormal principal directions in the laboratory frame.
% This function does not change the physical Hamiltonian or numerical step.
if nargin < 4, g_axes=eye(3); end
if isempty(coefficients_K), coefficients_K=zeros(0,1); end
validateattributes(coefficients_K,{'numeric'},{'real','finite'});
validateattributes(fields_T,{'numeric'},{'real','finite','2d','ncols',3,'nonempty'});
validateattributes(g,{'numeric'},{'real','finite','vector','numel',3});
validateattributes(g_axes,{'numeric'},{'real','finite','size',[3,3]});
if norm(g_axes*g_axes'-eye(3),'fro')>1e-10
    error('qspintn:EnergyScaleAxes','g principal directions must be orthonormal.');
end
u=QSpinThermalUnits();
effective=u.muB_over_kB_K_per_T*(fields_T*g_axes').*reshape(g,1,3);
% hypot avoids overflow in the norm when the finite inputs are large.
field_energies=hypot(hypot(effective(:,1),effective(:,2)),effective(:,3));
interaction_max=max([0;abs(coefficients_K(:))]);
field_max=max(field_energies);
reference=max(interaction_max,field_max);
zero_fallback=(reference==0);
if zero_fallback, reference=1; end
if ~isfinite(reference) || reference<=0
    error('qspintn:EnergyScale','Candidate energy scale must be positive and finite.');
end
scale=struct('schema_version','aqm.energy_scale.candidate_max.v1', ...
    'mode','candidate_max','energy_scale_K',reference, ...
    'interaction_max_K',interaction_max,'field_max_K',field_max, ...
    'field_energies_K',field_energies(:)', 'zero_fallback_used',zero_fallback, ...
    'zero_reference_K',1,'field_scope','all_evaluation_datasets');
end
