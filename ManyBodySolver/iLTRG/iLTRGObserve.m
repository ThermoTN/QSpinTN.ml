function [lambda, magnetization, diagnostics, env] = iLTRGObserve(Ta, La, Tb, Lb, Para, need_m, env)
%ILTRGOBSERVE Converged two-site bilayer transfer and spin expectation.
% lambda is the positive two-site transfer eigenvalue. magnetization is the
% bare spin projection along Para.Field.h, per spin; it contains no g factor.
% Para.transfer_tolerance and transfer_max_iterations default to 1e-12 and
% 10000. A failed convergence check is an error, never a usable measurement.
%
% The observable contractions are adapted from QSpinTN.ml/GetiLTRGM.m
% (Copyright (c) 2025 ThermoTN, MIT). Changes: strict fixed-point residuals,
% warm starts, imaginary-part checks, and explicit partition normalization.
% See the QSpinTN repository LICENSE for the retained MIT license.
if nargin < 6, need_m = true; end
if nargin < 7 || isempty(env), env = struct(); end
tol = setting(Para, 'transfer_tolerance', 1e-12);
max_iterations = setting(Para, 'transfer_max_iterations', 10000);
assert(isscalar(tol) && isfinite(tol) && tol > 0 && tol < 1, ...
    'qspintn:iltrg:InvalidTransferTolerance', 'Invalid transfer tolerance.');
assert(isscalar(max_iterations) && isfinite(max_iterations) && ...
    max_iterations >= 2 && max_iterations == floor(max_iterations), ...
    'qspintn:iltrg:InvalidTransferIterations', 'Invalid transfer iteration limit.');
left_seed = seed(env, 'left', size(La,1));
left_map = @(T) TransM(TransM(T, Ta, La, 'L2R'), Tb, Lb, 'L2R');
[Tl, lambda, left_info] = fixed_point(left_map, left_seed, tol, max_iterations, 'left');
env.left = Tl;
diagnostics = struct('tolerance',tol,'max_iterations',max_iterations, ...
    'left',left_info,'lambda_two_site',lambda,'converged',true);
magnetization = 0;
if ~need_m, return; end

right_seed = seed(env, 'right', size(La,1));
right_map = @(T) TransM(TransM(T, Tb, Lb, 'R2L'), Ta, La, 'R2L');
[Tr, right_lambda, right_info] = fixed_point(right_map, right_seed, tol, max_iterations, 'right');
env.right = Tr;
eigenvalue_difference = abs(log(right_lambda)-log(lambda));
if eigenvalue_difference > max(100*tol, 1e-10)
    error('qspintn:iltrg:TransferEigenvalueMismatch', ...
        'Left/right transfer log eigenvalues disagree by %.6g.', eigenvalue_difference);
end
diagnostics.right = right_info;
diagnostics.left_right_log_lambda_difference = eigenvalue_difference;

h = Para.Field.h(:).';
if numel(h) ~= 3 || ~all(isfinite(h)) || ~isreal(h) || norm(h) == 0
    error('qspintn:iltrg:InvalidObservableDirection', ...
        'Spin projection requires a real, nonzero three-component field.');
end
h = h/norm(h);
[Sx,Sy,Sz,spin] = physical_spin_operators(Para);
MOp = h(1)*Sx+h(2)*Sy+h(3)*Sz;

% This double layer represents rho(beta/2)*rho(beta/2)^dagger. Insert the
% observable on the evolved output physical leg, then normalize by the
% SAME complete AB transfer with no insertion. Using only <Tl,Tr> would
% leave an erroneous factor lambda unless the gauge happened to give 1.
Tap = contract(Ta,3,MOp,2,[1,2,4,3]);
Tmp = inserted_transfer(Tl,Tap,Ta,La);
Tmp = TransM(Tmp,Tb,Lb,'L2R');
Ma = contract(Tmp,[1,2],Tr,[1,2]);
Tbp = contract(Tb,3,MOp,2,[1,2,4,3]);
Tmp = TransM(Tl,Ta,La,'L2R');
Tmp = inserted_transfer(Tmp,Tbp,Tb,Lb);
Mb = contract(Tmp,[1,2],Tr,[1,2]);
full_transfer = left_map(Tl);
partition = contract(full_transfer,[1,2],Tr,[1,2]);
overlap = contract(Tl,[1,2],Tr,[1,2]);
partition = positive_real(partition, 'partition contraction', max(100*tol,1e-10));
normalization_difference = abs(partition-lambda*overlap)/partition;
if normalization_difference > max(100*tol,1e-10)
    error('qspintn:iltrg:TransferNormalizationMismatch', ...
        'Full transfer and eigenvalue partition normalizations differ by %.6g.', normalization_difference);
end
raw_m = (Ma+Mb)/(2*partition);
imaginary_m = abs(imag(raw_m));
if ~isfinite(raw_m) || imaginary_m > max(100*tol,1e-10)*max(1,abs(real(raw_m)))
    error('qspintn:iltrg:ComplexMagnetization', ...
        'Spin expectation is nonfinite or has a significant imaginary part (%.6g).', imaginary_m);
end
magnetization = real(raw_m);
if abs(magnetization) > spin+max(1000*tol,1e-9)
    error('qspintn:iltrg:UnphysicalMagnetization', ...
        'Spin expectation %.16g exceeds the spin bound %.16g.', magnetization,spin);
end
diagnostics.partition = partition;
diagnostics.partition_relative_difference = normalization_difference;
diagnostics.magnetization_imaginary = imaginary_m;
end

function [Sx,Sy,Sz,spin]=physical_spin_operators(Para)
sites=setting(Para,'physical_sites_per_tensor',1);
if isequal(sites,1)
    [Sx,Sy,Sz]=SpinOp(Para.d);
    spin=(Para.d-1)/2;
elseif isequal(sites,2) && Para.d==4
    % The d=4 local space is two physical spin-1/2 sites, NOT spin 3/2.
    % Insert their average; all resulting moments remain per physical spin.
    [sx,sy,sz]=SpinOp(2);
    Sx=(kron(sx,eye(2))+kron(eye(2),sx))/2;
    Sy=(kron(sy,eye(2))+kron(eye(2),sy))/2;
    Sz=(kron(sz,eye(2))+kron(eye(2),sz))/2;
    spin=.5;
else
    error('qspintn:iltrg:BlockProfile','Unsupported physical-spin observable representation.');
end
end

function T = inserted_transfer(T,inserted,original,L)
T = contract(T,1,L,1);
T = contract(T,1,conj(L),1);
T = contract(T,1,inserted,1);
T = contract(T,[1,3,4],conj(original),[1,3,4]);
end

function [T,lambda,info] = fixed_point(map,T,tol,max_iterations,direction)
previous_log_lambda = NaN;
stable_steps = 0;
maximum_hermiticity_error = 0;
for iteration = 1:max_iterations
    next = map(T);
    if any(~isfinite(next(:))) || norm(next,'fro') == 0
        error('qspintn:iltrg:InvalidTransfer', '%s transfer is zero or nonfinite.',direction);
    end
    next_norm = norm(next,'fro');
    hermiticity_error = norm(next-next','fro')/next_norm;
    maximum_hermiticity_error = max(maximum_hermiticity_error,hermiticity_error);
    if hermiticity_error > max(100*tol,1e-10)
        error('qspintn:iltrg:NonHermitianTransfer', ...
            '%s positive transfer lost Hermiticity: %.6g.',direction,hermiticity_error);
    end
    rayleigh = sum(conj(T(:)).*next(:))/sum(abs(T(:)).^2);
    lambda = positive_real(rayleigh,[direction,' eigenvalue'],max(100*tol,1e-10));
    residual = norm(next-lambda*T,'fro')/next_norm;
    log_lambda_change = abs(log(lambda)-previous_log_lambda);
    if residual <= tol && log_lambda_change <= tol
        stable_steps = stable_steps+1;
    else
        stable_steps = 0;
    end
    if stable_steps >= 2
        info = struct('converged',true,'iterations',iteration,'residual',residual, ...
            'log_lambda_change',log_lambda_change,'lambda',lambda, ...
            'maximum_hermiticity_error',maximum_hermiticity_error);
        return;
    end
    previous_log_lambda = log(lambda);
    % The map is completely positive. Remove roundoff-level anti-Hermitian
    % noise only after checking it; Frobenius normalization is positive.
    T = (next+next')/2;
    T = T/norm(T,'fro');
end
error('qspintn:iltrg:TransferNotConverged', ...
    '%s transfer did not converge in %d iterations (residual %.6g, dloglambda %.6g, tolerance %.6g).', ...
    direction,max_iterations,residual,log_lambda_change,tol);
end

function T = seed(env,name,D)
T = eye(D);
if isstruct(env) && isfield(env,name)
    candidate = env.(name);
    if isequal(size(candidate),[D,D]) && all(isfinite(candidate(:))) && norm(candidate,'fro') > 0
        % Full-rank identity admixture prevents a warm start from entirely
        % missing a sector after the SVD changes its virtual basis.
        candidate = (candidate+candidate')/2;
        T = .95*candidate/norm(candidate,'fro')+.05*eye(D)/sqrt(D);
    end
end
T = T/norm(T,'fro');
end

function result = positive_real(value,label,tol)
if ~isscalar(value) || ~isfinite(value) || real(value) <= 0 || abs(imag(value)) > tol*abs(real(value))
    error('qspintn:iltrg:InvalidPositiveTransferValue', ...
        '%s must be finite, positive and real (real %.6g, imaginary %.6g).',label,real(value),imag(value));
end
result = real(value);
end

function value = setting(Para,name,fallback)
value = fallback;
if isfield(Para,name), value = Para.(name); end
end
