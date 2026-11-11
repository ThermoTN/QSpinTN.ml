function [LnZ, M, raw] = iTEBD(Para, T)
%ITEBD Natural-unit wrapper around the shared iLTRG evolution core.
% Optional T is checked against the native grid; no interpolation is done.
% Explicit Para.terms_ab/terms_ba are natural-unit N-by-3 Trotter cell tables.
if isfield(Para,'terms_ab')
    ab=Para.terms_ab;
    if isfield(Para,'terms_ba'), ba=Para.terms_ba; else, ba=ab; end
else
    [ab,ba]=feval(Para.Trotter_Name,Para);
end
if nargin>=2
    expected=1./(2*Para.tau*(1:Para.N_max)');
    if numel(T)~=numel(expected) || any(abs(T(:)-expected)>1e-10*max(1,abs(expected)))
        error('qspintn:iltrg:TemperatureGrid','iTEBD expects its native temperature grid.');
    end
end
[Hab,Hba]=InitHam_trotter(ab,ba,Para);
if isscalar(Hab), Hab=Hab*eye(Para.d^2); end
if isscalar(Hba), Hba=Hba*eye(Para.d^2); end
% Retain only numerical fields: project labels, handles and filenames must
% not become part of a tensor state identity or control numerical behavior.
core=struct('d',Para.d,'D_max',Para.D_max,'N_max',Para.N_max, ...
    'tau',Para.tau,'energy_scale_K',1,'Field',struct('h',Para.Field.h));
for name={'transfer_tolerance','transfer_max_iterations'}
    if isfield(Para,name{1}), core.(name{1})=Para.(name{1}); end
end
if isscalar(Para.Field.h) && Para.Field.h==0
    % GetField uses a scalar magnitude; at zero field its hdir can be NaN.
    core.Field.h=[0,0,0];
elseif isfield(Para,'g_fec') && isscalar(Para.Field.h)
    direction=Para.Field.hdir(:).'; h=zeros(1,3);
    for k=1:3
        axis=Para.g_fec.dir{k}(:).';
        h=h+Para.Field.h*Para.g_fec.g(k)*dot(axis,direction)*axis;
    end
    core.Field.h=h;
end
options=struct();
if isfield(Para,'fileID'), options.log_fid=Para.fileID; end
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
old_path=path; restore_path=onCleanup(@() path(old_path)); %#ok<NASGU>
addpath(fullfile(root,'UtilityFunc'),fullfile(root,'svd_lapack_interface'), ...
    fullfile(root,'ManyBodySolver','iLTRG'),fullfile(root,'ManyBodySolver','TensorFunction'));
[raw,~,~]=iLTRGEvolve(core,Hab,Hba,options);
LnZ=raw.lnZ; M=raw.spin_projection;
end
