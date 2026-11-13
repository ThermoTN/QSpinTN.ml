function [result, state] = iLTRGForward(request, options)
%ILTRGFORWARD Physical-unit, explicit-model infinite-chain thermodynamics.
% [result,state] = iLTRGForward(request,options) accepts MATLAB structs.
% Interactions are in K, field in T, tau in K^-1. AB and BA bonds both
% follow the physical left-to-right orientation. Optional terms_nnn adds
% uniform i -> i+2 interactions via two-physical-spin blocks. See README.md.
% This entry point does not create project directories, JSON, CSV or logs.
% Copyright (c) 2025 ThermoTN. MIT License; see repository LICENSE.
if nargin<2, options=struct(); end
request=validate_request(request);
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
old_path=path; restore_path=onCleanup(@() path(old_path)); %#ok<NASGU>
for relative={'ManyBodySolver/iLTRG','ManyBodySolver/iLTRG/RG_func', ...
        'ManyBodySolver/TensorFunction','UtilityFunc','svd_lapack_interface'}
    addpath(fullfile(root,relative{1}),'-begin');
end
units=QSpinThermalUnits();
Para=struct('d',round(2*request.spin+1),'D_max',request.D_max, ...
    'N_max',request.N_max,'tau',request.tau*request.energy_scale_K, ...
    'energy_scale_K',request.energy_scale_K, ...
    'physical_sites_per_tensor',1,'physical_spin',request.spin, ...
    'transfer_tolerance',request.transfer_tolerance, ...
    'transfer_max_iterations',request.transfer_max_iterations,'request',request);
Para.Field.h=units.muB_over_kB_K_per_T*(request.g.*request.field_T)/request.energy_scale_K;
t_ab=term_table(request.terms_ab,request.energy_scale_K);
t_ba=term_table(request.terms_ba,request.energy_scale_K);
t_nnn=term_table(request.terms_nnn,request.energy_scale_K);
if isempty(t_nnn)
    [Hab,Hba]=InitHam_trotter(t_ab,t_ba,Para);
    profile='nearest_neighbor_spin_half'; maximum_neighbor_distance=1;
else
    Para.d=4; Para.physical_sites_per_tensor=2;
    [Hab,Hba]=InitHam_J1J2(t_ab,t_ba,t_nnn,Para);
    profile='nn_nnn_spin_half_two_site_blocks'; maximum_neighbor_distance=2;
end
if isscalar(Hab), Hab=Hab*eye(Para.d^2); end
if isscalar(Hba), Hba=Hba*eye(Para.d^2); end
[raw,state,info]=iLTRGEvolve(Para,Hab,Hba,options);
physical=ConvertThermalObservables(raw,request.field_T,request.g);
for name={'Cm_J_per_mol_K','M_muB_per_spin','Chi_emu_per_mol_Oe','Chi_qspintn_native'}
    raw.(name{1})=physical.(name{1});
end
if strcmp(request.observable,'Cm')
    values=raw.Cm_J_per_mol_K; column='Cm_J_per_mol_K'; unit='J/(mol K)';
else
    values=raw.Chi_emu_per_mol_Oe; column='Chi_emu_per_mol_Oe'; unit='emu/mol/Oe';
end
mask=isfinite(values) & raw.T_K>=request.T_min*(1-1e-12);
temperature=raw.T_K(mask); values=values(mask);
[temperature,order]=sort(temperature); values=values(order);
summary=rmfield(info,{'sources','science'});
summary.api_version='qspintn.iltrg.forward.v1';
summary.observable=request.observable; summary.unit=unit;
summary.model_profile=profile;
summary.maximum_neighbor_distance=maximum_neighbor_distance;
summary.tensor_local_dimension=Para.d;
summary.physical_sites_per_tensor=Para.physical_sites_per_tensor;
summary.physical_sites_per_unit_cell=2*Para.physical_sites_per_tensor;
summary.molar_basis=units.molar_basis; summary.points=numel(temperature);
summary.tau_K_inv=request.tau; summary.requested_T_min_K=request.T_min;
summary.singular_value_relative_cutoff=1e-8;
summary.interpolation='none; centered Cm on native internal beta grid; no extrapolation';
if ~isempty(temperature), summary.temperature_range_K=[min(temperature),max(temperature)]; end
if ~isempty(raw.beta)
    summary.maximum_D_a=max(raw.actual_D_a); summary.maximum_D_b=max(raw.actual_D_b);
end
if any(values < -1e-8)
    summary.warning='Negative values retained for diagnosis; no clipping or smoothing.';
end
science=struct('api_version',summary.api_version,'request',request, ...
    'sources',info.sources,'units',units,'identity',info.identity);
result=struct('status',info.status,'raw',raw,'summary',summary, ...
    'data',struct('T_K',temperature,'values',values,'column',column,'unit',unit), ...
    'hamiltonian',struct('Hab',Hab,'Hba',Hba,'Para',Para), ...
    'units',units,'science',science);
end

function req=validate_request(req)
defaults=struct('spin',.5,'energy_scale_K',1,'transfer_tolerance',1e-12, ...
    'transfer_max_iterations',10000,'terms_nnn',[]);
required={'terms_ab','field_T','g','observable','tau','D_max','N_max','T_min'};
allowed=[fieldnames(defaults);required(:);{'terms_ba'}];
if ~isstruct(req) || ~isscalar(req) || ~all(ismember(fieldnames(req),allowed))
    error('qspintn:iltrg:Request','Expected scientific request fields only.');
end
for name=fieldnames(defaults)'
    if ~isfield(req,name{1}), req.(name{1})=defaults.(name{1}); end
end
for name=required
    if ~isfield(req,name{1}), error('qspintn:iltrg:Request','Missing %s.',name{1}); end
end
if ~isfield(req,'terms_ba'), req.terms_ba=req.terms_ab; end
for name={'energy_scale_K','tau','T_min','transfer_tolerance'}
    validateattributes(req.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
for name={'D_max','N_max','transfer_max_iterations'}
    validateattributes(req.(name{1}),{'numeric'},{'scalar','real','finite','positive','integer'});
end
validateattributes(req.spin,{'numeric'},{'scalar','real','finite','positive'});
if req.spin~=.5
    error('qspintn:iltrg:Request','This validated physical forward profile requires spin=0.5.');
end
validateattributes(req.field_T,{'numeric'},{'vector','numel',3,'real','finite'});
validateattributes(req.g,{'numeric'},{'vector','numel',3,'real','finite','positive'});
req.field_T=reshape(req.field_T,1,3); req.g=reshape(req.g,1,3);
if isstring(req.observable) && isscalar(req.observable), req.observable=char(req.observable); end
if ~ischar(req.observable) || ~ismember(req.observable,{'Cm','Chi'})
    error('qspintn:iltrg:Request','observable must be Cm or Chi.');
end
if strcmp(req.observable,'Chi') && nnz(req.field_T)~=1
    error('qspintn:iltrg:Request','Chi requires one nonzero principal-axis field.');
end
if strcmp(req.observable,'Cm') && req.N_max<3
    error('qspintn:iltrg:Request','Cm needs at least three evolution steps.');
end
% Stable field order makes the checkpoint independent of struct construction.
req=orderfields(req);
end

function out=term_table(terms,scale)
if isempty(terms), out=cell(0,3); return; end
if ~isstruct(terms) || ~all(isfield(terms,{'operators','coefficient_K'})) || ...
        ~all(ismember(fieldnames(terms),{'operators','coefficient_K'}))
    error('qspintn:iltrg:Terms','Terms accept only operators and coefficient_K; their range is defined by the request field.');
end
out=cell(numel(terms),3);
for k=1:numel(terms)
    operators=terms(k).operators;
    if isstring(operators), operators=cellstr(operators); end
    if ~iscell(operators) || numel(operators)~=2 || ...
            ~all(ismember(operators,{'x','y','z','Sx','Sy','Sz'}))
        error('qspintn:iltrg:Terms','Each bilinear term needs two x/y/z operators.');
    end
    validateattributes(terms(k).coefficient_K,{'numeric'},{'scalar','real','finite'});
    for j=1:2
        op=char(operators{j}); if numel(op)==1, op=['S',op]; end
        out{k,j}=op;
    end
    out{k,3}=terms(k).coefficient_K/scale;
end
end
