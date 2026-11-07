function [state,energy,iterations] = DMRGLanczos(apply,state,Para)
% Shared three-term Lanczos, without reorthogonalization.
% Follow tanTRG_Spin/DMRG/LanczosGS: reuse the local effective operator,
% diagonalize the small tridiagonal matrix and stop on Ritz-energy change.
% Keep the existing two-site tolerance scaling, Para.tol/(10*Para.L).
shape = size(state); x = state(:); nx = norm(x);
assert(isfinite(nx) && nx>0,'DMRG:LanczosState','Invalid local initial state.');
K = min(Para.DK,numel(x)); q = cell(1,K); q{1} = x/nx;
tridiagonal = zeros(K); previous = Inf;
for k=1:K
    z = apply(reshape(q{k},shape)); z = z(:);
    diagonal = real(q{k}'*z);
    assert(isfinite(diagonal),'DMRG:LanczosEnergy','Nonfinite local energy.');
    tridiagonal(k,k) = diagonal;
    [vectors,values] = eig(tridiagonal(1:k,1:k),'vector');
    [energy,j] = min(real(values)); coefficients = vectors(:,j);
    scale = max(eps,abs(energy));
    % As in the original two-site solver/reference, stop also on a Ritz
    % energy increase (loss of numerical monotonicity in the Krylov space).
    if k>1 && (previous-energy)/scale*Para.L < Para.tol/10
        break
    end
    if k==K, break; end
    z = z-diagonal*q{k};
    if k>1, z = z-tridiagonal(k-1,k)*q{k-1}; end
    beta = norm(z);
    assert(isfinite(beta),'DMRG:LanczosResidual','Nonfinite Lanczos vector.');
    % Handle exact breakdown already at k=1, without dividing by zero.
    if beta <= Para.tol/10*scale, break; end
    tridiagonal(k,k+1) = beta; tridiagonal(k+1,k) = beta;
    q{k+1} = z/beta; previous = energy;
end
x = coefficients(1)*q{1};
for j=2:k, x = x+coefficients(j)*q{j}; end
nx = norm(x);
assert(isfinite(nx) && nx>0,'DMRG:LanczosState','Invalid Ritz state.');
state = reshape(x/nx,shape); iterations = k;
end
