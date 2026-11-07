function [T,filename] = LoadDMRGInitialMPS(filename)
% Load only the previous MPS; its Hamiltonian, parameters and environments
% do not belong to the new solve. Prefer the standard saved Rslt.T format.
T = [];
if isstring(filename)
    assert(isscalar(filename) && ~ismissing(filename),'DMRG:InitialMPSFile', ...
        'DMRGInitialMPSFile must be a scalar filename.');
    filename = char(filename);
end
assert(ischar(filename) && (isempty(filename) || isrow(filename)), ...
    'DMRG:InitialMPSFile','DMRGInitialMPSFile must be a filename or an empty string.');
if isempty(filename), return; end
assert(isfile(filename),'DMRG:InitialMPSFile','Initial MPS file does not exist: %s',filename);
[ok,attributes] = fileattrib(filename);
assert(ok,'DMRG:InitialMPSFile','Cannot resolve initial MPS file: %s',filename);
filename = attributes.Name;
variables = whos('-file',filename); names = {variables.name};
if ismember('Rslt',names)
    saved = load(filename,'Rslt');
    if isstruct(saved.Rslt) && isscalar(saved.Rslt) && isfield(saved.Rslt,'T')
        T = saved.Rslt.T;
    end
end
if isempty(T) && ismember('T',names)
    saved = load(filename,'T'); T = saved.T;
end
assert(~isempty(T),'DMRG:InitialMPSFile', ...
    'File must contain Rslt.T or T. Generate the previous result with saveMPS=true: %s',filename);
end
