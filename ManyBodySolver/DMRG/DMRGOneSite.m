function [T,energy,EE,info] = DMRGOneSite(Para,H,env,handoff)
% Finite, open-MPS single-site DMRG; H may describe periodic physical bonds.
% A(left,right,physical), W(left,right,out,in), environment(bra,ket,MPO).
L = Para.L; d = Para.d;
sweepTol = GetDMRGSweepTol(Para);
validateattributes(L,{'numeric'},{'scalar','integer','>=',2});
validateattributes(d,{'numeric'},{'scalar','integer','>=',2});
validateattributes(Para.MCrit,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(Para.tol,{'numeric'},{'scalar','positive','finite'});
validateattributes(Para.DK,{'numeric'},{'scalar','integer','>=',2,'finite'});
validateattributes(Para.DMRGStepMax,{'numeric'},{'scalar','positive','integer','finite'});
assert(numel(H.A)==L,'MPO length does not match Para.L.');
W = H.A;
W{1} = reshape(W{1},1,[],d,d); W{L} = reshape(W{L},[],1,d,d);
sweepOffset = 0;
if nargin<4 || isempty(handoff)
    A = initialState(L,d,Para.MCrit);
else
    % A completed backward two-site sweep has its center at site 1. All
    % other tensors are already right isometries; reshape only the ends.
    A = handoff.T.A;
    A{1} = reshape(A{1},1,[],d); A{L} = reshape(A{L},[],1,d);
    sweepOffset = handoff.Sweep;
end
right = 1;
for i=L:-1:2
    right = extendRight(right,A{i},W{i}); env.put(i-1,'R',right);
end
initialEnergy = centerEnergy(A{1},W{1},right,H.lgnorm);
if sweepOffset>0
    assert(abs(initialEnergy-handoff.Energy)<=1e-9*max(1,abs(handoff.Energy)), ...
        'DMRG:HandoffEnergy','Rebuilt energy %.16g differs from two-site energy %.16g.', ...
        initialEnergy,handoff.Energy);
end
info = struct('InitialEnergy',initialEnergy,'Sweeps',0,'Converged',false, ...
    'SweepSeconds',[],'LocalIterations',[]);
bondEE = zeros(L-1,1); oldEE = 0; oldEnergy = NaN;
for sweep=1:Para.DMRGStepMax
    ticSweep = tic;
    localIterations = 0;
    left = 1;
    for i=1:L-1
        right = env.get(i,'R');
        [center,k] = optimize(A{i},left,W{i},right,Para);
        localIterations = localIterations+k;
        [A{i},A{i+1}] = shiftRight(center,A{i+1});
        left = extendLeft(left,A{i},W{i}); env.put(i,'L',left);
    end
    right = 1;
    for i=L:-1:2
        left = env.get(i-1,'L');
        [center,k] = optimize(A{i},left,W{i},right,Para);
        localIterations = localIterations+k;
        [A{i},A{i-1},bondEE(i-1)] = shiftLeft(center,A{i-1});
        right = extendRight(right,A{i},W{i}); env.put(i-1,'R',right);
    end
    energy = centerEnergy(A{1},W{1},right,H.lgnorm);
    assert(isfinite(energy),'DMRG:Energy','Nonfinite single-site energy at sweep %d.',sweepOffset+sweep);
    EE = max(bondEE);
    [converged,deltaS] = CheckDMRGSweep(energy,oldEnergy,EE,oldEE,sweepTol);
    fprintf('%d, En: %.16f, DeltaS: %g\n     Trunc. Er.: 0\n     EE: %.16f\n', ...
        sweepOffset+sweep,energy,deltaS,EE);
    info.Sweeps = sweep;
    info.SweepSeconds(sweep) = toc(ticSweep);
    info.LocalIterations(sweep) = localIterations;
    fprintf('Elapsed time is %.6f seconds.\n',info.SweepSeconds(sweep));
    if sweep>1 && converged
        info.Converged = true; break
    end
    oldEnergy = energy; oldEE = EE;
end
% Restore the two-index end tensors used by GetDMRGObs and existing readers.
A{1} = reshape(A{1},[],d); A{L} = reshape(A{L},[],d);
T = struct('A',{A},'lgnorm',0);
end

function A = initialState(L,d,maxBond)
% Single-site sweeps cannot create new bond sectors from a product state.
% Start at the requested rank, capped by the physical Hilbert space of each cut.
bonds = ones(1,L+1);
for b=2:L, bonds(b)=min(maxBond,d*bonds(b-1)); end
bonds = min(bonds,fliplr(bonds));
A = cell(L,1);
for i=1:L
    [Q,~] = qr(randn(bonds(i+1)*d,bonds(i)),0);
    A{i} = reshape(Q',bonds(i),bonds(i+1),d);
end
end

function energy = centerEnergy(center,W,right,logScale)
applied = action(1,center,W,right);
energy = real(center(:)'*applied(:))*exp(logScale);
end

function B = leftAction(left,A,W)
% Partial H_eff*A, ordered (bra-left, ket-right, MPO-right, physical-out).
B = tensorprod(left,A,2,1,NumDimensionsA=3);
B = tensorprod(B,W,[2,4],[1,4],NumDimensionsA=4);
end

function Y = action(left,A,W,right)
B = leftAction(left,A,W);
Y = permute(tensorprod(B,right,[2,3],[2,3],NumDimensionsA=4),[1,3,2]);
end

function E = extendLeft(left,A,W)
E = tensorprod(conj(A),leftAction(left,A,W),[1,3],[1,4],NumDimensionsA=3);
end

function E = extendRight(right,A,W)
E = extendLeft(right,permute(A,[2,1,3]),permute(W,[2,1,3,4]));
end

function [center,iterations] = optimize(center,left,W,right,Para)
% Like the reference DMRG_OS_update, fuse the local MPO with the left
% environment once. Cache matrix layouts outside the Krylov loop as well.
% Only this site's fused block is resident, also with the disk backend.
[dl,dr,d] = size(center); wr = size(W,2);
fused = tensorprod(left,W,3,1,NumDimensionsA=3);
fused = reshape(permute(fused,[1,4,3,2,5]),dl*d*wr,dl*d);
rightMatrix = reshape(permute(right,[2,3,1]),dr*wr,dr);
apply = @(x) fusedAction(x,fused,rightMatrix,dl,dr,d,wr);
[center,~,iterations] = DMRGLanczos(apply,center,Para);
end

function Y = fusedAction(A,fused,rightMatrix,dl,dr,d,wr)
B = fused*reshape(permute(A,[1,3,2]),dl*d,dr);
B = reshape(permute(reshape(B,dl*d,wr,dr),[1,3,2]),dl*d,dr*wr);
Y = permute(reshape(B*rightMatrix,dl,d,dr),[1,3,2]);
end

function [site,next] = shiftRight(center,next)
[dl,dr,d] = size(center);
matrix = reshape(permute(center,[1,3,2]),dl*d,dr);
[Q,R] = qr(matrix,0); newBond = size(Q,2);
site = permute(reshape(Q,dl,d,newBond),[1,3,2]);
next = reshape(R*reshape(next,dr,[]),newBond,size(next,2),d);
end

function [site,prev,entropy] = shiftLeft(center,prev)
[dl,dr,d] = size(center);
[U,S,V] = svd(reshape(center,dl,dr*d),'econ'); newBond = size(V,2);
site = reshape(V',newBond,dr,d);
prevLeft = size(prev,1);
mat = reshape(permute(prev,[1,3,2]),prevLeft*d,dl)*(U*S);
prev = permute(reshape(mat,prevLeft,d,newBond),[1,3,2]);
p = diag(S).^2; p = p(p>0)/sum(p);
entropy = -sum(p.*log(p));
end
