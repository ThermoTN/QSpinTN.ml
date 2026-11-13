function sources = iLTRGSourceIdentity()
%ILTRGSOURCEIDENTITY Relocatable content identity of the numerical profile.
% Record the resolved LAPACK fallback binary as well as its MATLAB interface.
% Absolute source paths are deliberately excluded from checkpoint identity.
root=canonical(fileparts(fileparts(fileparts(mfilename('fullpath')))));
core='ManyBodySolver/iLTRG'; tensor='ManyBodySolver/TensorFunction';
rg=[core,'/RG_func']; utility='UtilityFunc'; lapack_dir='svd_lapack_interface';
profile={'iLTRGForward',core; 'iLTRGEvolve',core; 'iLTRGObserve',core; ...
    'iLTRGSourceIdentity',core; 'InitHam_trotter',core; 'InitHam_J1J2',core; ...
    'InitId',rg; 'EvoGate',rg; 'EvoRG',rg; 'TransM',rg; ...
    'SpinOp',tensor; 'contract',tensor; 'svdT',tensor; ...
    'svd_lapack',lapack_dir; 'QSpinThermalUnits',utility; ...
    'ConvertThermalObservables',utility};
sources=struct('name',{},'path',{},'sha256',{});
for k=1:size(profile,1)
    name=profile{k,1}; expected=canonical(fullfile(root,profile{k,2},[name,'.m']));
    resolved=which(name);
    if isempty(resolved) || ~strcmp(canonical(resolved),expected)
        error('qspintn:iltrg:Profile','%s must resolve to %s.',name,expected);
    end
    sources(end+1)=record(name,expected,root); %#ok<AGROW>
end
expected=fullfile(root,lapack_dir,['lapack.',mexext]);
if ~isfile(expected), expected=fullfile(root,lapack_dir,'lapack.m'); end
resolved=which('lapack');
if isempty(resolved) || ~strcmp(canonical(resolved),canonical(expected))
    error('qspintn:iltrg:Profile','LAPACK fallback must resolve to %s.',expected);
end
sources(end+1)=record('lapack',canonical(expected),root);
% The .m fallback is also distributed when a platform binary is present.
matlab_fallback=fullfile(root,lapack_dir,'lapack.m');
if ~strcmp(canonical(expected),canonical(matlab_fallback))
    sources(end+1)=record('lapack_matlab_fallback',canonical(matlab_fallback),root);
end
end

function item=record(name,location,root)
fid=fopen(location,'rb');
if fid<0, error('qspintn:iltrg:Profile','Cannot read numerical source %s.',location); end
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
bytes=fread(fid,inf,'*uint8');
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(typecast(uint8(bytes(:)),'int8'));
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
item=struct('name',name,'path',location(numel(root)+2:end),'sha256',digest);
end
function value=canonical(value)
value=char(java.io.File(char(value)).getCanonicalPath());
end
