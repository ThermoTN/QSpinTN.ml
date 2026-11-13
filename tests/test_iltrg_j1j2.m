function report = test_iltrg_j1j2(runtime_root)
%TEST_ILTRG_J1J2 Independent checks of the blocked spin-1/2 J1-J2 API.
% Artifacts belong outside the source checkout. The references below use
% explicit Pauli matrices, a classical transfer matrix and finite-ring ED;
% none use the production block-Hamiltonian constructor or its spin operators.
% These modest-D/high-T checks do not certify arbitrary low-T convergence.
package=fileparts(fileparts(mfilename('fullpath')));
if nargin<1, runtime_root=tempdir; end
assert(java.io.File(runtime_root).isAbsolute(),'Use an absolute runtime directory.');
runtime_root=char(java.io.File(runtime_root).getCanonicalPath());
assert(~strcmp(runtime_root,package) && ~startsWith(runtime_root,[package,filesep]), ...
    'Artifacts must be outside the source checkout.');
if ~isfolder(runtime_root), mkdir(runtime_root); end
root=tempname(runtime_root); mkdir(root);
old_path=path; restore_path=onCleanup(@() path(old_path)); %#ok<NASGU>
addpath(fullfile(package,'ManyBodySolver','iLTRG'));
assert(strcmp(which('iLTRGForward'),fullfile(package,'ManyBodySolver','iLTRG','iLTRGForward.m')));
log_fid=fopen(fullfile(root,'evolution.log'),'w'); assert(log_fid>0);
close_log=onCleanup(@() fclose(log_fid)); %#ok<NASGU>
options=struct('log_fid',log_fid);
R=8.31446; muB_over_kB=9.274e-24/(0.08617*1.602e-22);
emu_per_muB=6.022e23*9.274e-24*1e3;
zero=term('z','z',0);
base=struct('spin',.5,'terms_ab',zero,'terms_ba',zero,'terms_nnn',zero, ...
    'energy_scale_K',1,'field_T',[0,0,0],'g',[1.5,1.8,2], ...
    'observable','Cm','tau',.01,'D_max',8,'N_max',24,'T_min',2);

% Bond range is selected exclusively by terms_ab/ba/nnn. Silently ignoring
% metadata such as distance=3 would simulate a different Hamiltonian.
invalid=base; invalid.terms_nnn.distance=3;
validation_checks.nnn_distance=rejection(invalid,options,'qspintn:iltrg:Terms');
invalid=base; invalid.terms_nnn.unknown_option=true;
validation_checks.nnn_unknown_field=rejection(invalid,options,'qspintn:iltrg:Terms');
invalid=base; invalid.terms_ab.distance=2;
validation_checks.nn_ab_extra_field=rejection(invalid,options,'qspintn:iltrg:Terms');
invalid=base; invalid.terms_ba.unknown_option=true;
validation_checks.nn_ba_extra_field=rejection(invalid,options,'qspintn:iltrg:Terms');
invalid=base; invalid.spin=1;
validation_checks.unsupported_spin=rejection(invalid,options,'qspintn:iltrg:Request');
invalid=base; invalid.terms_nnn.operators={'x','raising'};
validation_checks.invalid_operator=rejection(invalid,options,'qspintn:iltrg:Terms');
invalid=base; invalid.terms_nnn.coefficient_K=1+1i;
validation_checks.complex_coefficient=rejection(invalid,options,'MATLAB:expectedReal');

% A three-block periodic ring has six physical spins. Comparing complete
% matrices catches missing intra-block half weights, duplicate fields, the
% AB/BA distinction and directed nearest/next-nearest DM orientation.
request=base;
request.terms_ab=[term('x','x',.8),term('y','y',.7),term('z','z',.6), ...
    term('x','y',.19),term('y','x',-.19),term('x','z',.13)];
request.terms_ba=[term('x','x',.51),term('y','y',.42),term('z','z',.31), ...
    term('x','y',-.27),term('y','x',.27),term('z','x',-.17)];
request.terms_nnn=[term('x','x',.33),term('z','z',.29), ...
    term('x','y',.21),term('y','x',-.21),term('x','z',.11)];
request.field_T=[.15,-.11,.23]; request.energy_scale_K=2.7;
request.N_max=3; request.tau=.002;
output=iLTRGForward(request,options);
numerical_sources=output.science.sources;
assert(output.hamiltonian.Para.d==4);
assert(output.hamiltonian.Para.physical_sites_per_tensor==2);
assert(isequal(size(output.hamiltonian.Hab),[16,16]));
assert(norm(output.hamiltonian.Hab-output.hamiltonian.Hba,'fro')<1e-14);
direct=physical_ring(request,6,muB_over_kB);
blocked=zeros(64);
for site=1:3
    blocked=blocked+embed_pair(output.hamiltonian.Hab,4,3,site,mod(site,3)+1);
end
blocked=blocked*request.energy_scale_K;
matrix_error=norm(blocked-direct,'fro');
assert(matrix_error<3e-13,'Blocked periodic Hamiltonian disagrees with physical bonds: %.3g',matrix_error);
matrix_check=struct('physical_sites',6,'blocks',3,'frobenius_error_K',matrix_error, ...
    'includes_alternating_NN',true,'includes_directed_NNN_DM',true,'includes_oblique_field',true);

% Explicit zero NNN coefficients force the grouped representation, testing
% per-physical-spin lnZ, Cm and moment normalization on all principal axes.
free_checks=cell(1,3);
for axis=1:3
    request=base; request.field_T(axis)=.2; request.observable='Chi';
    request.tau=.02; request.N_max=24; request.T_min=1;
    output=iLTRGForward(request,options); beta=output.raw.beta(:);
    h=muB_over_kB*request.g(axis)*request.field_T(axis);
    exact_log=log(2*cosh(beta*h/2)); exact_m=.5*tanh(beta*h/2);
    log_error=max(abs(output.raw.lnZ(:)-exact_log));
    spin_error=max(abs(output.raw.spin_projection(:)-exact_m));
    chi_error=max(abs(output.raw.Chi_emu_per_mol_Oe(:)- ...
        emu_per_muB*request.g(axis)*exact_m/(1e4*request.field_T(axis))));
    c_error=heat_difference_error(output.raw,exact_log);
    assert(log_error<2e-8 && spin_error<2e-8 && chi_error<3e-7);
    assert(c_error<2e-6,'Blocked free-spin Cm normalization mismatch.');
    assert(max(abs(output.raw.M_muB_per_spin(:)))<=request.g(axis)/2+1e-10);
    assert(max(abs(output.raw.beta(:)-2*request.tau*(1:request.N_max)'))<1e-14);
    free_checks{axis}=struct('axis',axis,'lnZ_per_spin_max_error',log_error, ...
        'spin_max_error',spin_error,'Chi_emu_per_mol_Oe_max_error',chi_error, ...
        'C_natural_centered_difference_error',c_error);
end

% Commuting J1z-J2z is solved by an independent 4x4 pair-state transfer
% matrix at every native beta. The field and both exchanges are nonzero.
classical=base; classical.terms_ab=term('z','z',.83);
classical.terms_ba=classical.terms_ab; classical.terms_nnn=term('z','z',.27);
classical.field_T=[0,0,.17]; classical.observable='Chi';
classical.tau=.025; classical.N_max=32; classical.T_min=.5;
classical_output=iLTRGForward(classical,options);
h=muB_over_kB*classical.g(3)*classical.field_T(3);
[exact_log,exact_m]=classical_reference(classical_output.raw.beta,.83,.27,h);
classical_log_error=max(abs(classical_output.raw.lnZ(:)-exact_log));
classical_spin_error=max(abs(classical_output.raw.spin_projection(:)-exact_m));
classical_c_error=heat_difference_error(classical_output.raw,exact_log);
assert(classical_log_error<1e-7 && classical_spin_error<1e-7, ...
    'Classical J1-J2 infinite-chain transfer reference mismatch.');
assert(classical_c_error<5e-6,'Classical J1-J2 centered Cm mismatch: %.3g',classical_c_error);
assert(max(abs(classical_output.raw.Chi_qspintn_native-4*pi*classical_output.raw.Chi_emu_per_mol_Oe))<1e-11);
classical_check=struct('lnZ_per_spin_max_error',classical_log_error, ...
    'spin_max_error',classical_spin_error,'C_natural_centered_difference_error',classical_c_error, ...
    'reference','Exact infinite-chain pair-state transfer matrix, not finite ED');

scaled=classical; scaled.energy_scale_K=3.7;
scaled_output=iLTRGForward(scaled,options);
scale_error=raw_difference(scaled_output.raw,classical_output.raw);
assert(scale_error<2e-7,'Changing E0 at fixed physical tau changed observables: %.3g',scale_error);

% With J1=0, even and odd spins are two independent quantum chains. Their
% per-physical-spin thermodynamics equals that of one unblocked NN chain.
split=base; split.terms_nnn=[term('x','x',.7),term('y','y',.7),term('z','z',.9)];
split.field_T=[0,0,.1]; split.observable='Chi'; split.D_max=24;
split.N_max=16; split.T_min=3;
split_output=iLTRGForward(split,options);
single=rmfield(split,'terms_nnn'); single.terms_ab=split.terms_nnn;
single.terms_ba=single.terms_ab; single.D_max=16;
single_output=iLTRGForward(single,options);
split_log_error=max(abs(split_output.raw.lnZ-single_output.raw.lnZ));
split_spin_error=max(abs(split_output.raw.spin_projection-single_output.raw.spin_projection));
split_c_error=finite_difference(split_output.raw.C_natural,single_output.raw.C_natural);
assert(split_log_error<2e-5 && split_spin_error<2e-5 && split_c_error<5e-4, ...
    'Decoupled J2 chains do not reproduce a single NN chain within the stated numerical tolerance.');
split_check=struct('lnZ_max_difference',split_log_error,'spin_max_difference',split_spin_error, ...
    'C_natural_max_difference',split_c_error,'blocked_D',split.D_max,'single_chain_D',single.D_max, ...
    'scope','High-temperature finite-D comparison; truncations differ between grouped and ungrouped representations');

% Noncommuting uniform Heisenberg J1-J2: compare two independently assembled
% periodic ED sizes, then change D and tau separately at common native betas.
quantum=base; quantum.terms_ab=heisenberg(1); quantum.terms_ba=quantum.terms_ab;
quantum.terms_nnn=heisenberg(.4); quantum.field_T=[0,0,.2];
quantum.observable='Chi'; quantum.tau=.01; quantum.N_max=16; quantum.T_min=3;
quantum.D_max=16; coarse=iLTRGForward(quantum,options);
quantum.D_max=24; larger=iLTRGForward(quantum,options);
refined=quantum; refined.tau=.005; refined.N_max=32;
fine=iLTRGForward(refined,options);
beta=fine.raw.beta(:);
[log8,m8,c8]=quantum_reference(quantum,8,beta,muB_over_kB);
[log10,m10,c10]=quantum_reference(quantum,10,beta,muB_over_kB);
finite_size=struct('lnZ_max_difference',max(abs(log8-log10)), ...
    'spin_max_difference',max(abs(m8-m10)),'C_natural_max_difference',max(abs(c8-c10)));
assert(finite_size.lnZ_max_difference<2e-6 && finite_size.spin_max_difference<2e-6 ...
    && finite_size.C_natural_max_difference<2e-5,'ED references are not size-stable in the chosen high-T interval.');
ed_log_error=max(abs(fine.raw.lnZ-log10));
ed_spin_error=max(abs(fine.raw.spin_projection-m10));
ed_c_error=heat_difference_error(fine.raw,log10);
assert(ed_log_error<5e-5 && ed_spin_error<5e-5 && ed_c_error<1e-3, ...
    'High-temperature quantum J1-J2 disagrees with independent ED beyond the declared tolerance.');
d_changes=struct('lnZ_max_difference',max(abs(coarse.raw.lnZ-larger.raw.lnZ)), ...
    'spin_max_difference',max(abs(coarse.raw.spin_projection-larger.raw.spin_projection)), ...
    'C_natural_max_difference',finite_difference(coarse.raw.C_natural,larger.raw.C_natural));
common=2:2:refined.N_max;
tau_changes=struct('lnZ_max_difference',max(abs(larger.raw.lnZ-fine.raw.lnZ(common))), ...
    'spin_max_difference',max(abs(larger.raw.spin_projection-fine.raw.spin_projection(common))), ...
    'C_natural_max_difference',finite_difference(larger.raw.C_natural,fine.raw.C_natural(common)));
assert(d_changes.lnZ_max_difference<5e-5 && d_changes.spin_max_difference<5e-5 ...
    && d_changes.C_natural_max_difference<1e-3);
assert(tau_changes.lnZ_max_difference<5e-5 && tau_changes.spin_max_difference<5e-5 ...
    && tau_changes.C_natural_max_difference<1e-3);
quantum_check=struct('J1_K',1,'J2_K',.4,'beta_range_K_inv',[min(beta),max(beta)], ...
    'ED_sizes',[8,10],'ED_finite_size',finite_size,'lnZ_max_error_vs_ED10',ed_log_error, ...
    'spin_max_error_vs_ED10',ed_spin_error,'C_natural_centered_error_vs_ED10',ed_c_error, ...
    'Cm_J_per_mol_K_centered_error_vs_ED10',R*ed_c_error, ...
    'D_comparison',d_changes,'tau_comparison',tau_changes, ...
    'settings',[struct('D',16,'tau_K_inv',.01),struct('D',24,'tau_K_inv',.01),struct('D',24,'tau_K_inv',.005)], ...
    'scope','High-temperature numerical agreement; not a low-temperature or global convergence certification');

% Grouped NN-only is a distinct Trotter partition of the same Hamiltonian.
nn=quantum; nn.terms_nnn=zero;
blocked_nn=iLTRGForward(nn,options); nn=rmfield(nn,'terms_nnn');
unblocked_nn=iLTRGForward(nn,options);
nn_error=struct('lnZ_max_difference',max(abs(blocked_nn.raw.lnZ-unblocked_nn.raw.lnZ)), ...
    'spin_max_difference',max(abs(blocked_nn.raw.spin_projection-unblocked_nn.raw.spin_projection)), ...
    'C_natural_max_difference',finite_difference(blocked_nn.raw.C_natural,unblocked_nn.raw.C_natural));
assert(nn_error.lnZ_max_difference<5e-5 && nn_error.spin_max_difference<5e-5 && nn_error.C_natural_max_difference<1e-3);

% Durable and in-memory continuations must reproduce the same grouped
% trajectory. A changed J2 coefficient is a new scientific checkpoint.
restart_options=options; restart_options.checkpoint_path=fullfile(root,'restart.mat');
restart_options.checkpoint_interval=3; restart_options.pause_after_steps=5;
[partial,state]=iLTRGForward(classical,restart_options);
assert(strcmp(partial.status,'paused') && state.iteration==5);
restart_options=rmfield(restart_options,'pause_after_steps');
resumed=iLTRGForward(classical,restart_options);
assert(resumed.summary.resumed_from_iteration==5 && strcmp(resumed.status,'completed'));
memory_options=options; memory_options.state=state;
memory=iLTRGForward(classical,memory_options);
restart_error=raw_difference(resumed.raw,classical_output.raw);
memory_error=raw_difference(memory.raw,classical_output.raw);
assert(restart_error<1e-10 && memory_error<1e-10);
changed=classical; changed.terms_nnn.coefficient_K=.28;
rejected=false;
try, iLTRGForward(changed,restart_options);
catch err, rejected=strcmp(err.identifier,'qspintn:iltrg:CheckpointMismatch'); end
assert(rejected,'Changed J2 must not resume the old tensor checkpoint.');

report=struct('status','passed','artifact_dir',root,'block_matrix',matrix_check, ...
    'request_validation',validation_checks,'numerical_sources',numerical_sources, ...
    'free_spins',{free_checks},'classical_J1z_J2z',classical_check, ...
    'energy_scale_raw_max_difference',scale_error,'decoupled_J2_chains',split_check, ...
    'quantum_J1_J2',quantum_check,'blocked_NN_limit',nn_error, ...
    'file_resume_max_difference',restart_error,'memory_resume_max_difference',memory_error, ...
    'changed_J2_checkpoint_rejected',rejected, ...
    'normalization','All thermodynamic quantities are per physical spin, not per d=4 block', ...
    'scope','Independent matrix/transfer/finite-ED checks and modest-D high-T convergence comparisons only');
fid=fopen(fullfile(root,'j1j2_report.json'),'w'); assert(fid>0);
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');
fprintf('QSPINTN_ILTRG_J1J2_TEST_PASS %s\n',fullfile(root,'j1j2_report.json'));
end

function item=rejection(request,options,expected_identifier)
actual_identifier='';
try
    iLTRGForward(request,options);
catch err
    actual_identifier=err.identifier;
end
assert(strcmp(actual_identifier,expected_identifier), ...
    'Invalid request returned "%s" instead of "%s".',actual_identifier,expected_identifier);
item=struct('rejected',true,'error_identifier',actual_identifier);
end
function item=term(a,b,c)
item=struct('operators',{{a,b}},'coefficient_K',c);
end
function terms=heisenberg(j)
terms=[term('x','x',j),term('y','y',j),term('z','z',j)];
end
function S=pauli()
S={[0,1;1,0]/2,[0,1i;-1i,0]/2,diag([-1,1])/2};
end
function index=axis_index(name)
if numel(name)==2 && name(1)=='S', name=name(2); end
index=find('xyz'==name,1); assert(~isempty(index));
end
function [H,M]=physical_ring(request,n,muB_over_kB)
S=pauli(); identity=speye(2); operators=cell(n,3);
for i=1:n
    for axis=1:3
        value=1;
        for j=1:n
            if i==j, value=kron(value,sparse(S{axis})); else, value=kron(value,identity); end
        end
        operators{i,axis}=value;
    end
end
H=sparse(2^n,2^n); M=sparse(2^n,2^n);
h=muB_over_kB*(request.g.*request.field_T);
for i=1:n
    if mod(i,2)==1, nearest=request.terms_ab; else, nearest=request.terms_ba; end
    for distance=1:2
        if distance==1, terms=nearest; else, terms=request.terms_nnn; end
        j=mod(i+distance-1,n)+1;
        for k=1:numel(terms)
            a=axis_index(terms(k).operators{1}); b=axis_index(terms(k).operators{2});
            H=H+terms(k).coefficient_K*(operators{i,a}*operators{j,b});
        end
    end
    for axis=1:3, H=H-h(axis)*operators{i,axis}; end
    M=M+operators{i,3};
end
H=full(H); M=full(M);
assert(norm(H-H','fro')<1e-12);
end
function embedded=embed_pair(local,d,n,i,j)
dimension=d^n; states=zeros(dimension,n);
for site=1:n, states(:,site)=mod(floor((0:dimension-1)'/d^(n-site)),d)+1; end
other=setdiff(1:n,[i,j]); embedded=zeros(dimension);
for row=1:dimension
    for col=1:dimension
        if all(states(row,other)==states(col,other))
            embedded(row,col)=local((states(row,i)-1)*d+states(row,j), ...
                                     (states(col,i)-1)*d+states(col,j));
        end
    end
end
end
function [logZ,m]=classical_reference(beta,j1,j2,h)
states=[-.5,-.5;-.5,.5;.5,-.5;.5,.5]; logZ=zeros(numel(beta),1); m=logZ;
for k=1:numel(beta)
    transfer=zeros(4); insertion=zeros(4);
    for a=1:4
        for b=1:4
            if states(a,2)==states(b,1)
                transfer(a,b)=exp(-beta(k)*(j1*states(a,1)*states(a,2) ...
                    +j2*states(a,1)*states(b,2)-h*states(a,1)));
                insertion(a,b)=transfer(a,b)*states(a,1);
            end
        end
    end
    [v,e]=eig(transfer); [lambda,index]=max(real(diag(e))); right=real(v(:,index));
    [v,e]=eig(transfer.'); [~,index]=min(abs(diag(e)-lambda)); left=real(v(:,index));
    logZ(k)=log(lambda); m(k)=(left.'*insertion*right)/(lambda*(left.'*right));
end
end
function [logZ,m,c]=quantum_reference(request,n,beta,muB_over_kB)
[H,M]=physical_ring(request,n,muB_over_kB); [vectors,energy]=eig(H,'vector');
energy=real(energy); shift=min(energy);
weights=exp(-beta(:)*(energy(:)-shift).'); partition=sum(weights,2); weights=weights./partition;
logZ=(log(partition)-beta(:)*shift)/n;
spin=real(sum(conj(vectors).*(M*vectors),1)).'; m=weights*spin/n;
mean_energy=weights*energy;
c=beta(:).^2.*(weights*(energy.^2)-mean_energy.^2)/n;
end
function error_value=heat_difference_error(raw,exact_log)
beta=raw.beta(:); mid=2:numel(beta)-1; step=beta(2)-beta(1);
exact=beta(mid).^2.*(exact_log(mid-1)+exact_log(mid+1)-2*exact_log(mid))/step^2;
error_value=max(abs(raw.C_natural(mid)-exact));
end
function error_value=finite_difference(a,b)
assert(isequal(size(a),size(b))); mask=isfinite(a)&isfinite(b);
assert(any(mask(:))); error_value=max(abs(a(mask)-b(mask)));
end
function error_value=raw_difference(a,b)
error_value=0;
for field={'beta','T_K','lnZ','spin_projection','C_natural','M_muB_per_spin', ...
        'Chi_emu_per_mol_Oe','Chi_qspintn_native','actual_D_a','actual_D_b'}
    x=a.(field{1}); y=b.(field{1});
    assert(isequal(size(x),size(y)) && isequal(isnan(x),isnan(y)));
    mask=isfinite(x)&isfinite(y);
    if any(mask(:)), error_value=max(error_value,max(abs(x(mask)-y(mask)))); end
end
end
