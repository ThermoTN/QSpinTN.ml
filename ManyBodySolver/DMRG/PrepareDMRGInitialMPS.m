function T = PrepareDMRGInitialMPS(T,Para)
% Validate and right-canonicalize a supplied MPS without SVD truncation.
% Package layout: endpoints (bond,physical), bulk (left,right,physical).
validateattributes(Para.L,{'numeric'},{'scalar','integer','>=',2,'finite'});
validateattributes(Para.d,{'numeric'},{'scalar','integer','>=',2,'finite'});
validateattributes(Para.MCrit,{'numeric'},{'scalar','integer','positive','finite'});
assert(isstruct(T) && isscalar(T) && isfield(T,'A') && iscell(T.A) && ...
    numel(T.A)==Para.L,'DMRG:InitialMPSShape','Initial MPS must contain Para.L tensors in T.A.');
L = Para.L; d = Para.d; A = T.A(:); previousBond = 1;
for site=1:L
    tensor = A{site};
    assert((isa(tensor,'double') || isa(tensor,'single')) && ~issparse(tensor) && ...
        ~isempty(tensor) && all(isfinite(tensor),'all'), ...
        'DMRG:InitialMPSShape','Initial MPS tensor %d must be a finite dense floating-point array.',site);
    if site==1
        assert(ismatrix(tensor) && size(tensor,2)==d,'DMRG:InitialMPSShape', ...
            'First MPS tensor must have shape [bond, Para.d].');
        A{site} = reshape(tensor,1,size(tensor,1),d);
    elseif site==L
        assert(ismatrix(tensor) && size(tensor,2)==d,'DMRG:InitialMPSShape', ...
            'Last MPS tensor must have shape [bond, Para.d].');
        A{site} = reshape(tensor,size(tensor,1),1,d);
    else
        assert(ndims(tensor)==3 && size(tensor,3)==d,'DMRG:InitialMPSShape', ...
            'MPS tensor %d must have shape [left, right, Para.d].',site);
    end
    assert(size(A{site},1)==previousBond,'DMRG:InitialMPSShape', ...
        'MPS bond dimensions do not match at site %d.',site);
    previousBond = size(A{site},2);
    assert(previousBond<=Para.MCrit,'DMRG:InitialMPSBond', ...
        'Initial MPS bond at site %d is %d, above MCrit=%d. Increase MCrit or supply a smaller MPS.', ...
        site,previousBond,Para.MCrit);
end
% Right QR moves the center to site 1 and preserves the represented state.
for site=L:-1:2
    dl = size(A{site},1); dr = size(A{site},2);
    [Q,R] = qr(reshape(A{site},dl,dr*d)',0);
    rank = size(Q,2); A{site} = reshape(Q',rank,dr,d);
    left = size(A{site-1},1);
    matrix = reshape(permute(A{site-1},[1,3,2]),left*d,dl)*R';
    A{site-1} = permute(reshape(matrix,left,d,rank),[1,3,2]);
end
scale = norm(A{1}(:));
assert(isfinite(scale) && scale>0,'DMRG:InitialMPSNorm','Initial MPS has zero or nonfinite norm.');
A{1} = A{1}/scale;
A{1} = reshape(A{1},[],d); A{L} = reshape(A{L},[],d);
% Global amplitude is irrelevant after normalization; discard old metadata.
T = struct('A',{A},'lgnorm',0);
end
