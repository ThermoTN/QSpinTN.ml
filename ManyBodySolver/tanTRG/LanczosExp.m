function [B, Info] = LanczosExp(F_handle, A, delta, varargin)
% Usage: [B, Info]= LanczosExp(F_handle, initState A, delta [, opts])
% calculate exp(H*delta)*A by Lanczos exponentiation
% F_handle is the function handle, must linearly hermitide act on A
%   'K': dimension of Krylov space (default = 32)
%   'tol': tolerance of truncating Krylov space (default = 1e-8)
%   'plus': function handle Y =  A + B, the same below
%   'mtimes':  Y = a*A (a is a scalar)
%   'norm':  y = |A|
%   'innerprod': y = <A,B>
%   'sum': y = sum(...), accelerate the summation of the last step
% warning: Class of A must have inner product space structure, i.e. method
% 'plus', 'minus', 'mtimes' (linear combination) and 'norm', 'innerpord'
% must be given if they has not been overloaded for class of A (default)
% note that inner product is implemented by polarization identity, which is
% very expensive, we prefer to give this parameter to accelerate

P = inputParser;
P.addParameter('K', 32, @(x) isnumeric(x) && x >= 1);
P.addParameter('tol', 1e-8, @(x) isnumeric(x) && x > 0);
P.addParameter('plus', @(x,y) plus(x,y), @(x) isa(x, 'function_handle'));
P.addParameter('mtimes', @(x,y) mtimes(x,y), @(x) isa(x, 'function_handle'));
P.addParameter('norm', @(x) norm(x), @(x) isa(x, 'function_handle'));
P.addParameter('innerprod', [], @(x) isa(x, 'function_handle'));
P.addParameter('sum', [], @(x) isa(x, 'function_handle'));
P.parse(varargin{:});
name_fields = fields(P.Results);
for i = 1:numel(fields(P.Results))
    eval(sprintf('tmp = P.Results.%s;', name_fields{i}));
    if isa(tmp, 'function_handle')
        name = sprintf('F_%s',  name_fields{i});
    else
        name =  name_fields{i};
    end
    eval([name,'=tmp;']);
end

Info.name = 'LanczosExp';
Info.opts.K = K;
Info.opts.tol = tol;

% use plus and mtimes to get minus
F_minus = @(x,y) F_plus(x, F_mtimes(-1, y));
% use norm to get innerprod (polarization identity of complex inner product space)
if ~exist('F_innerprod','var')
    F_innerprod = @(x,y) (F_norm(F_plus(x, y))^2 - F_norm(F_minus(x, y))^2 ...
        + 1i*F_norm(F_plus(x, F_mtimes(1i, y)))^2 ...
        - 1i*F_norm(F_minus(x, F_mtimes(1i, y)))^2)/4;
end

% Lanczos basis
q = cell(1, K+1);
% tridiagonal form of the linear operator on Lanczos basis
T = zeros(K+1);
% initial state
q{1} = A;
% normalize
Norm = F_norm(q{1});
q{1} = F_mtimes(Norm^(-1), q{1});
for k = 1:K
    % H*q_k
    q{k+1} = F_handle(q{k});
    % <q_k|H|q_k>
    T(k,k) = F_innerprod(q{k},q{k+1});
    % main step of Lanczos, update next Lanczos base
    % q_{k+1} = q_{k+1} - T(k,k)*q_k - T(k-1,k)*q_{k-1}
    q{k+1} = F_minus(q{k+1}, F_mtimes(T(k, k), q{k}));
    if k > 1
        q{k+1} = F_minus(q{k+1}, F_mtimes(T(k-1,k),q{k-1}));
    end
    % update off-diagonal
    T(k, k + 1) = F_norm(q{k + 1});
    T(k + 1, k) = T(k, k + 1);
    
    % calculate coefficient of expm
    coef = expm(delta*T(1:k, 1:k))*[1; zeros(k - 1, 1)];

    % check linear independence
    if abs(T(k, k + 1)/T(1, 1)) < tol
        Info.Err = abs(T(k, k + 1)/T(1, 1));
        break;
    end
    
    % check if more dimension if needed
    if abs(coef(end)/norm(coef)) < tol
        Info.Err = abs(coef(end)/norm(coef));
        break;
    end

    % normalize
    q{k + 1} = F_mtimes(T(k,k+1)^(-1),q{k + 1});  
end

% use the last coefficient to estimate error
if k == K
    Info.Err = abs(coef(end));
end

% update Lanczos info
Info.K = k; % actual Krylov dimension
Info.T = T(1:k, 1:k); % the tridiagonal matrix
Info.Ek = eig(Info.T, 'vector');
Info.coef = coef;
% fprintf('K: %d\n', Info.K)

% calculate target states
if ~exist('F_sum','var')
    % use plus
    B = F_mtimes(coef(1), q{1});
    for j = 2:Info.K
        B = F_plus(B, F_mtimes(coef(j), q{j}));
    end
else
    % use sum to accelerate (VariSum of global MPS for example)
    S_list = cell(1, Info.K);
    for j = 1:Info.K
        S_list{j} = F_mtimes(coef(j), q{j});
    end
    B = F_sum(S_list{:});
end

% remember the norm factor
B = F_mtimes(Norm, B);

end
