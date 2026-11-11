function report = test_iltrg_forward(runtime_root)
%TEST_ILTRG_FORWARD Independent physical checks of the package iLTRG API.
%   report = test_iltrg_forward(absolute_runtime_directory)
% Run from outside the source checkout; isolate MATLAB_PREFDIR and TMPDIR in
% the calling process. No agent, JSON request, or autoQMagen code is required.
% Small bond dimensions here do not establish production low-T convergence.
package = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1, runtime_root = tempdir; end
assert(java.io.File(runtime_root).isAbsolute(), 'An absolute runtime directory is required.');
runtime_root = char(java.io.File(runtime_root).getCanonicalPath());
assert(~strcmp(runtime_root,package) && ~startsWith(runtime_root,[package,filesep]), ...
    'Test artifacts must be outside the QSpinTN source checkout.');
if ~isfolder(runtime_root), mkdir(runtime_root); end
root = tempname(runtime_root); mkdir(root);
old_path = path; path_cleanup = onCleanup(@() path(old_path)); %#ok<NASGU>
addpath(fullfile(package,'ManyBodySolver','iLTRG'), ...
    fullfile(package,'ManyBodySolver','iLTRG','RG_func'), ...
    fullfile(package,'ManyBodySolver','TensorFunction'), ...
    fullfile(package,'UtilityFunc'),fullfile(package,'svd_lapack_interface'));
assert(strcmp(which('iLTRGForward'),fullfile(package,'ManyBodySolver','iLTRG','iLTRGForward.m')));

% Independent ascending-Sz basis, matching QSpinTN's basis convention.
S = {[0,1;1,0]/2, [0,1i;-1i,0]/2, diag([-1,1])/2};
H = kron(S{1},S{1}) + kron(S{2},S{2}) + 1.5*kron(S{3},S{3}) ...
    + .5*(kron(S{1},S{2})-kron(S{2},S{1}));
expected_spectrum = sort([.375;.375;-.375-sqrt(1.25)/2;-.375+sqrt(1.25)/2]);
assert(norm(H-H','fro')<1e-14 && any(imag(H(:))~=0));
assert(max(abs(sort(real(eig(H)))-expected_spectrum))<1e-14);

% Isolate all four right-operator branches of InitHam_trotter. Spin 1
% distinguishes Sz^2 from the identity; this validates only the local
% constructor and does not extend the validated spin-1/2 forward profile.
spin_one_x = [0,1,0;1,0,1;0,1,0]/sqrt(2);
spin_one_y = [0,1i,0;-1i,0,1i;0,-1i,0]/sqrt(2);
spin_one_z2 = diag([1,0,1]);
constructor_para = struct('d',3,'Field',struct('h',[0,0,0]));
constructor_checks = cell(1,4); constructor_index = 0;
for operator_index=1:2
    right_names = {'Sz2','Id'}; left_names = {'Sx','Sy'};
    left_matrices = {spin_one_x,spin_one_y};
    right_matrices = {spin_one_z2,eye(3)};
    terms_local = {left_names{operator_index},right_names{operator_index},.37};
    zero_local = {'Sz','Sz',0};
    expected_local = .37*kron(left_matrices{operator_index},right_matrices{operator_index});
    for bond_index=1:2
        if bond_index==1
            [actual_ab,actual_ba] = InitHam_trotter(terms_local,zero_local,constructor_para);
            constructor_error = max(norm(actual_ab-expected_local,'fro'),norm(actual_ba,'fro'));
            bond_name = 'AB';
        else
            [actual_ab,actual_ba] = InitHam_trotter(zero_local,terms_local,constructor_para);
            constructor_error = max(norm(actual_ab,'fro'),norm(actual_ba-expected_local,'fro'));
            bond_name = 'BA';
        end
        assert(constructor_error<1e-14, ...
            'Right %s on %s must change the right operator only.',right_names{operator_index},bond_name);
        assert(norm(actual_ab-actual_ab','fro')<1e-14 && norm(actual_ba-actual_ba','fro')<1e-14);
        constructor_index = constructor_index+1;
        constructor_checks{constructor_index} = struct('bond',bond_name, ...
            'right_operator',right_names{operator_index},'spin',1, ...
            'local_H_max_error',constructor_error,'scope','local Hamiltonian construction only');
    end
end

kb = 8.617e-2*1.602e-22; muB = 9.274e-24;
R = 8.31446; emu_per_muB = 5584.8028;
% The common converter must keep native susceptibility separate from CGS.
unit_raw = struct('C_natural',[1;2],'spin_projection',[.1;.25]);
physical = ConvertThermalObservables(unit_raw,[0,0,.1],[1.5,1.5,2]);
expected_unit_chi = emu_per_muB*(2*unit_raw.spin_projection)/1000;
assert(max(abs(physical.Cm_J_per_mol_K-R*unit_raw.C_natural))<1e-13);
assert(max(abs(physical.M_muB_per_spin-2*unit_raw.spin_projection))<1e-13);
assert(max(abs(physical.Chi_emu_per_mol_Oe-expected_unit_chi))<1e-12);
assert(max(abs(physical.Chi_qspintn_native-4*pi*expected_unit_chi))<1e-11);
negative_field = ConvertThermalObservables(unit_raw,[0,0,-.1],[1.5,1.5,2]);
assert(max(abs(negative_field.M_muB_per_spin+physical.M_muB_per_spin))<1e-13);
assert(max(abs(negative_field.Chi_emu_per_mol_Oe-physical.Chi_emu_per_mol_Oe))<1e-12);
zero_field = ConvertThermalObservables(unit_raw,[0,0,0],[1.5,1.5,2]);
assert(all(isnan(zero_field.Chi_emu_per_mol_Oe)) && all(zero_field.M_muB_per_spin==0));
unit_checks = struct('native_to_cgs_factor',1/(4*pi), ...
    'CGS_and_native_distinct',true,'negative_field_consistent',true, ...
    'zero_field_M_over_H_undefined',true);

terms = [term('x','x',1),term('y','y',1),term('z','z',1.5), ...
    term('x','y',.5),term('y','x',-.5)];
zero_terms = term('x','x',0);
base = struct('spin',.5,'energy_scale_K',1,'terms_ab',terms,'terms_ba',terms, ...
    'field_T',[0,0,0],'g',[1.5,1.5,2],'observable','Cm', ...
    'tau',.01,'D_max',16,'N_max',53,'T_min',1);

% Noninteracting spins are exact at any Trotter step. The same half-field
% on both AB and BA bonds must give a full single-site Zeeman energy once.
free_checks = cell(1,2);
for k=1:2
    axis = [1,3]; axis = axis(k);
    request = base; request.terms_ab = zero_terms; request.terms_ba = zero_terms;
    request.field_T(axis) = .1; request.observable = 'Chi'; request.D_max = 8;
    tag = sprintf('free_axis_%d',axis);
    [summary,raw,output] = invoke(request,root,tag);
    assert(strcmp(summary.status,'completed'));
    beta = raw.beta(:); h = muB/kb*request.g(axis)*.1;
    exact_logZ = log(2*cosh(beta*h/2));
    exact_spin = .5*tanh(beta*h/2);
    exact_moment = request.g(axis)*exact_spin;
    exact_chi = emu_per_muB*exact_moment/(1e4*.1);
    spin_error = max(abs(raw.spin_projection(:)-exact_spin));
    moment_error = max(abs(raw.M_muB_per_spin(:)-exact_moment));
    chi_error = max(abs(raw.Chi_emu_per_mol_Oe(:)-exact_chi));
    assert(spin_error<2e-8,'Free-spin expectation mismatch on axis %d: %.3g',axis,spin_error);
    assert(moment_error<4e-8 && chi_error<3e-7,'Experimental magnetization/Chi conversion mismatch.');
    native_para = struct('Model',struct('ES',request.energy_scale_K));
    native_units = GetUnitCon(native_para);
    natural_h = h/request.energy_scale_K;
    native_chi = request.g(axis)^2 * raw.spin_projection(:)/natural_h * native_units.Chi_con;
    native_to_cgs_error = max(abs(native_chi/(4*pi)-raw.Chi_emu_per_mol_Oe(:)));
    assert(native_to_cgs_error<3e-12, 'Native QSpinTN chi and CGS chi must differ by 4*pi.');
    assert(strcmp(output.data.unit,'emu/mol/Oe'),'CGS data must carry an explicit unit.');
    expected_field_H = -.5*h*(kron(S{axis},eye(2))+kron(eye(2),S{axis}));
    ham = output.hamiltonian;
    field_error = max(norm(ham.Hab-expected_field_H,'fro'),norm(ham.Hba-expected_field_H,'fro'));
    assert(field_error<1e-13,'g/field must be applied once and split evenly across the two bonds.');
    logZ_check = check_logZ(raw,exact_logZ);
    exact_C = (beta*h/2).^2 ./ cosh(beta*h/2).^2;
    heat_check = check_heat(raw,exact_logZ,exact_C,R,2e-6);
    assert_export(output,request.T_min,request.N_max,raw.Chi_emu_per_mol_Oe,'Chi_emu_per_mol_Oe');
    free_checks{k} = struct('axis',axis,'field_energy_K',h,'spin_max_error',spin_error, ...
        'moment_muB_per_spin_max_error',moment_error,'Chi_emu_per_mol_Oe_max_error',chi_error, ...
        'field_H_max_error_K',field_error,'native_to_cgs_max_error',native_to_cgs_error, ...
        'logZ',logZ_check,'heat_capacity',heat_check);
end

% Observable normalization must be independent of a tensor's scalar gauge.
% This directly detects a missing transfer-eigenvalue denominator, including
% the imaginary Sy physical leg, without an SVD or evolution approximation.
gauge_checks = cell(1,3);
gauge_beta = .7; gauge_h = .3;
gauge_old_path = path;
gauge_path_cleanup = onCleanup(@() path(gauge_old_path));
addpath(fullfile(package,'ManyBodySolver','TensorFunction'), ...
    fullfile(package,'ManyBodySolver','iLTRG','RG_func'));
for axis=1:3
    rho_half = expm(gauge_beta*gauge_h*S{axis}/2);
    tensor = reshape(rho_half,[1,1,2,2]);
    h = zeros(1,3); h(axis) = gauge_h;
    para = struct('d',2,'Field',struct('h',h));
    [lambda,moment] = iLTRGObserve(tensor,1,tensor,1,para,true);
    [scaled_lambda,scaled_moment] = iLTRGObserve(2*tensor,1,tensor,1,para,true);
    expected_lambda = (2*cosh(gauge_beta*gauge_h/2))^2;
    expected_moment = .5*tanh(gauge_beta*gauge_h/2);
    assert(abs(lambda-expected_lambda)<1e-12 && abs(scaled_lambda-4*lambda)<1e-12);
    moment_error = max(abs([moment,scaled_moment]-expected_moment));
    assert(moment_error<1e-12,'Magnetization must be invariant under tensor scalar gauge.');
    gauge_checks{axis} = struct('axis',axis,'lambda',lambda, ...
        'scaled_lambda',scaled_lambda,'spin_max_error',moment_error);
end
clear gauge_path_cleanup;

% The existing package entry also uses the common evolution. A pure Chi
% request must return a defined LnZ on the same un-interpolated temperature
% grid, and the natural M must agree with the analytic free-spin limit.
legacy = struct('d',2,'D_max',8,'N_max',base.N_max,'tau',base.tau);
legacy.terms_ab = {'Sx','Sx',0}; legacy.terms_ba = legacy.terms_ab;
legacy.Field.h = [0,0,muB/kb*2*.1];
legacy_result = GetiLTRGRslt(legacy,'Chi');
assert(numel(legacy_result.T)==numel(legacy_result.beta));
assert(numel(legacy_result.T)==numel(legacy_result.LnZ));
assert(numel(legacy_result.T)==numel(legacy_result.M));
assert(numel(legacy_result.T)==numel(legacy_result.Cm));
assert(all(isfinite(legacy_result.T)) && all(isfinite(legacy_result.M)));
assert(max(abs(legacy_result.T(:).*legacy_result.beta(:)-1))<1e-12);
legacy_exact_M = .5*tanh(legacy_result.beta(:)*legacy.Field.h(3)/2);
legacy_moment_error = max(abs(legacy_result.M(:)-legacy_exact_M));
assert(legacy_moment_error<2e-8);
legacy_lnZ = check_logZ(struct('lnZ',legacy_result.LnZ), ...
    log(2*cosh(legacy_result.beta(:)*legacy.Field.h(3)/2)));
native_beta = 2*legacy.tau*(1:legacy.N_max);
for k=1:numel(legacy_result.beta)
    assert(min(abs(native_beta-legacy_result.beta(k)))<1e-12, ...
        'The package result must not add interpolated or extrapolated temperatures.');
end
% GetField's historical zero-field representation may have a scalar h and
% an undefined direction. Both forms must mean the same zero Hamiltonian.
legacy_zero = legacy;
legacy_zero.Field = struct('h',0,'hdir',[NaN,NaN,NaN]);
scalar_zero = GetiLTRGRslt(legacy_zero,'Cm');
legacy_zero.Field = struct('h',[0,0,0]);
vector_zero = GetiLTRGRslt(legacy_zero,'Cm');
zero_field_error = 0;
for name={'T','beta','LnZ','M','Cm'}
    a=scalar_zero.(name{1}); b=vector_zero.(name{1});
    assert(isequal(size(a),size(b)) && all(isfinite(a)) && all(isfinite(b)));
    zero_field_error = max(zero_field_error,max(abs(a(:)-b(:))));
end
assert(zero_field_error<1e-12,'Scalar and vector zero fields must produce the same result.');
assert(max(abs(scalar_zero.Cm))<1e-8 && max(abs(scalar_zero.M))<1e-12);
assert(max(abs(scalar_zero.LnZ-log(2)))<1e-8, ...
    'Noninteracting zero-field spin-1/2 must have ln Z per spin = ln 2.');
legacy_checks = struct('pure_Chi_defined',true,'aligned_native_grid',true, ...
    'spin_max_error',legacy_moment_error,'logZ',legacy_lnZ, ...
    'scalar_vector_zero_field_max_error',zero_field_error);


% One AB dimer and zero BA coupling: exact independent four-level dimers.
% This checks the complex DM gate and per-spin normalization without any
% finite-chain or thermodynamic-limit approximation in the reference.
dimer = base; dimer.terms_ba = zero_terms; dimer.tau = .005;
dimer.N_max = 103; dimer.D_max = 8;
[dimer_summary,dimer_raw,dimer_dir] = invoke(dimer,root,'independent_dimers');
assert(strcmp(dimer_summary.status,'completed'));
ham = dimer_dir.hamiltonian;
ham_error = max(norm(ham.Hab-H,'fro'),norm(ham.Hba,'fro'));
assert(ham_error<1e-13,'Directed DM interaction map has the wrong local Hamiltonian.');
assert(norm(ham.Hab-ham.Hab','fro')<1e-13 && any(imag(ham.Hab(:))~=0));
assert(max(abs(sort(real(eig(ham.Hab)))-expected_spectrum))<1e-13);
beta = dimer_raw.beta(:);
shift = min(expected_spectrum);
weights = exp(-beta*(expected_spectrum-shift)');
partition = sum(weights,2); weights = weights./partition;
exact_logZ = .5*(log(partition)-beta*shift);
mean_E = weights*expected_spectrum;
exact_C = .5*beta.^2.*(weights*(expected_spectrum.^2)-mean_E.^2);
dimer_logZ = check_logZ(dimer_raw,exact_logZ);
dimer_heat = check_heat(dimer_raw,exact_logZ,exact_C,R,2e-4);
assert_export(dimer_dir,dimer.T_min,dimer.N_max,R*dimer_raw.C_natural,'Cm_J_per_mol_K');
assert(strcmp(dimer_dir.data.unit,'J/(mol K)'));
scaled_dimer = dimer; scaled_dimer.energy_scale_K = 2.75;
[scaled_summary,scaled_raw,scaled_result] = invoke(scaled_dimer,root,'scaled_dimers');
assert(strcmp(scaled_summary.status,'completed'));
assert(max(abs(scaled_raw.beta-dimer_raw.beta))<1e-14);
assert(norm(scaled_result.hamiltonian.Hab*scaled_dimer.energy_scale_K-H,'fro')<1e-13);
scaled_logZ = check_logZ(scaled_raw,exact_logZ);
scaled_heat = check_heat(scaled_raw,exact_logZ,exact_C,R,2e-4);
assert_export(scaled_result,dimer.T_min,dimer.N_max,R*scaled_raw.C_natural,'Cm_J_per_mol_K');

% Actual uniform XXZ+Dz model, at short beta and a small bond dimension.
% Its purpose is finite output, endpoint and directed-bond consistency;
% D=16 here is not evidence of production convergence at low temperature.
uniform = base; uniform.observable = 'Chi'; uniform.field_T = [.1,0,0];
uniform.N_max = 28; uniform.T_min = 2;
[uniform_summary,uniform_raw,uniform_dir] = invoke(uniform,root,'uniform_complete');
assert(strcmp(uniform_summary.status,'completed'));
ham = uniform_dir.hamiltonian;
hx = muB/kb*uniform.g(1)*.1;
expected_uniform_H = H-.5*hx*(kron(S{1},eye(2))+kron(eye(2),S{1}));
assert(norm(ham.Hab-expected_uniform_H,'fro')<1e-13);
assert(norm(ham.Hba-expected_uniform_H,'fro')<1e-13, ...
    'BA is a left-to-right physical bond too: do not reverse its DM sign.');
assert(all(isfinite(uniform_raw.lnZ(:))) && all(isfinite(uniform_raw.spin_projection(:))));
assert(isreal(uniform_raw.lnZ) && isreal(uniform_raw.spin_projection));
assert(all(abs(uniform_raw.spin_projection(:))<=.5+1e-10),'Spin-1/2 magnetization exceeded saturation.');
assert(all(isfinite(uniform_raw.C_natural(2:end-1))));
assert_export(uniform_dir,uniform.T_min,uniform.N_max,uniform_raw.Chi_emu_per_mol_Oe,'Chi_emu_per_mol_Oe');

% A controlled partial invocation must save the exact tensor state; a
% second invocation continues at step 10, and a third does no evolution.
[paused,partial_raw,~,partial_state] = invoke(uniform,root,'uniform_resumed',struct('pause_after_steps',9));
assert(strcmp(paused.status,'paused') && numel(partial_raw.beta)==9);
[resumed,resumed_raw] = invoke(uniform,root,'uniform_resumed');
assert(strcmp(resumed.status,'completed') && resumed.resumed_from_iteration==9);
resume_error = compare_raw(resumed_raw,uniform_raw);
assert(resume_error<2e-12,'Checkpoint resume changed the tensor evolution: %.3g',resume_error);
[memory_result,~] = iLTRGForward(uniform,struct('state',partial_state));
assert(strcmp(memory_result.status,'completed') && memory_result.summary.resumed_from_iteration==9);
memory_resume_error = compare_raw(memory_result.raw,uniform_raw);
assert(memory_resume_error<2e-12,'Memory-state resume changed the tensor evolution.');
[reused,reused_raw] = invoke(uniform,root,'uniform_resumed');
assert(strcmp(reused.status,'completed') && reused.resumed_from_iteration==uniform.N_max);
assert(isequaln(reused_raw,resumed_raw),'Completed checkpoint reuse must not recompute raw results.');
changed = uniform; changed.terms_ab(1).coefficient_K = 1.01;
rejected = false;
try
    invoke(changed,root,'uniform_resumed');
catch exception
    rejected = strcmp(exception.identifier,'qspintn:iltrg:CheckpointMismatch');
    if ~rejected, rethrow(exception); end
end
assert(rejected,'A checkpoint for different physical parameters must be rejected.');

report = struct('passed',true,'run_root',root,'matlab_version',version, ...
    'constants',struct('kB_J_per_K',kb,'muB_J_per_T',muB,'R_J_per_mol_K',R, ...
        'emu_per_mol_per_muB_per_spin',emu_per_muB), ...
    'directed_two_spin_spectrum_K',expected_spectrum,'local_H_max_error_K',ham_error, ...
    'right_operator_constructor_checks',{constructor_checks}, ...
    'free_spin_checks',{free_checks},'observable_gauge_checks',{gauge_checks}, ...
    'unit_contract',unit_checks,'legacy_entry',legacy_checks, ...
    'dimer_logZ',dimer_logZ,'dimer_heat_capacity',dimer_heat, ...
    'scaled_energy_logZ',scaled_logZ,'scaled_energy_heat_capacity',scaled_heat, ...
    'uniform_steps',uniform.N_max,'uniform_minimum_exported_T_K',uniform.T_min, ...
    'resumed_from_iteration',resumed.resumed_from_iteration,'resume_raw_max_error',resume_error, ...
    'memory_resume_raw_max_error',memory_resume_error, ...
    'completed_checkpoint_reused',true,'changed_parameters_rejected',rejected, ...
    'scope',['Small-D package API validation against free spins and independent dimers; ', ...
        'uniform-chain finiteness and checkpoint validation, not a low-temperature convergence study.']);
write_json(fullfile(root,'iltrg_report.json'),report);
fprintf('QSPINTN_ILTRG_TEST_PASS %s\n',fullfile(root,'iltrg_report.json'));
end

function item = term(left,right,coefficient)
item = struct('operators',{{left,right}},'coefficient_K',coefficient);
end

function [summary,raw,result,state] = invoke(request,root,tag,options)
if nargin < 4, options = struct(); end
options.checkpoint_path = fullfile(root,[tag,'_checkpoint.mat']);
options.checkpoint_interval = 7;
[result,state] = iLTRGForward(request,options);
summary = result.summary; summary.status = result.status;
raw = result.raw;
end

function check = check_logZ(raw,expected)
delta = raw.lnZ(:)-expected(:);
offset = mean(delta);
shape_error = max(abs(delta-offset));
assert(shape_error<5e-8,'log Z temperature dependence mismatch: %.3g',shape_error);
% An additive normalization constant cancels in energy/Cm derivatives.
% Record it explicitly so it cannot be confused with an absolute entropy.
check = struct('additive_offset_per_spin',offset,'max_error_after_constant_offset',shape_error, ...
    'absolute_normalization_agrees',abs(offset)<5e-8);
end

function check = check_heat(raw,exact_logZ,exact_C,R,physical_tolerance)
beta = raw.beta(:); step = diff(beta);
assert(max(abs(step-step(1)))<1e-12 && all(beta>0));
assert(max(abs(raw.T_K(:)-1./beta))<1e-12);
middle = 2:(numel(beta)-1);
exact_difference = beta(middle).^2.*(exact_logZ(middle-1)+exact_logZ(middle+1) ...
    -2*exact_logZ(middle))/step(1)^2;
actual = raw.C_natural(:);
difference_error = max(abs(actual(middle)-exact_difference));
physical_error = max(abs(R*actual(middle)-R*exact_C(middle)));
assert(difference_error<2e-6,'Cm must agree with the independent centered logZ derivative: %.3g',difference_error);
assert(physical_error<physical_tolerance,'Analytical Cm mismatch: %.3g J/(mol spin K)',physical_error);
check = struct('C_natural_max_error_vs_exact_centered_difference',difference_error, ...
    'Cm_J_per_mol_spin_K_max_error_vs_exact_derivative',physical_error, ...
    'analytic_tolerance_J_per_mol_spin_K',physical_tolerance);
end

function assert_export(result,T_min,N_max,expected,value_column)
data = result.data;
assert(strcmp(data.column,value_column));
assert(~isempty(data.T_K) && numel(data.T_K)<=N_max);
values = data.values;
assert(all(isfinite(data.T_K)) && all(isfinite(values)));
assert(all(diff(data.T_K)>0),'Exported physical temperatures must increase.');
assert(abs(min(data.T_K)-T_min)<1e-12,'The requested minimum T must be present without extrapolation.');
for k=1:numel(data.T_K)
    [distance,index] = min(abs(result.raw.T_K(:)-data.T_K(k)));
    assert(distance<1e-11,'Exported values must remain on native temperature points.');
    assert(abs(values(k)-expected(index))<2e-10,'Experimental-unit result conversion mismatch.');
end
end

function error_value = compare_raw(actual,expected)
fields = {'beta','T_K','lnZ','spin_projection','C_natural','M_muB_per_spin', ...
    'Chi_emu_per_mol_Oe','actual_D_a','actual_D_b'};
error_value = 0;
for k=1:numel(fields)
    a = actual.(fields{k}); b = expected.(fields{k});
    assert(isequal(size(a),size(b)) && isequal(isnan(a),isnan(b)));
    finite = isfinite(a) & isfinite(b);
    assert(all(finite(:)|isnan(a(:))));
    if any(finite(:)), error_value = max(error_value,max(abs(a(finite)-b(finite)))); end
end
end

function write_json(path,value)
fid = fopen(path,'w'); assert(fid~=-1,'Cannot write test artifact: %s',path);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid,jsonencode(value,'PrettyPrint',true),'char');
end
