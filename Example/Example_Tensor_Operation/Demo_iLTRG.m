clear all; clc
addpath(genpath('../../ManyBodySolver'))

% Set parameter 
Para.tau = 0.01;
Para.D = 200;
Para.d = 2;
Para.maxsetp = 5/(Para.tau * 2);

% Get H
[Sx, Sy, Sz, Id] = SpinOp(Para.d);
H = kron(Sx, Sx) + kron(Sy, Sy);

% Get U
U = expm(-Para.tau * H);
U = reshape(U, [Para.d, Para.d, Para.d, Para.d]);

% Init T
La = 1;
Lb = 1;
Ta = reshape(Id, [1, 1, Para.d, Para.d]);
Tb = reshape(Id, [1, 1, Para.d, Para.d]);

lgnorm = 0;
for it = 1:Para.maxsetp
    beta(it) = Para.tau * 2 * it;
    fprintf('it: %6d| beta: %6f| ', it, beta(it));
    [ Ta, Lb, Tb, Ns ] = RG( La, Ta, Lb, Tb, U, Para);
    lgnorm = lgnorm + log(Ns);
    [ Tb, La, Ta, Ns ] = RG( Lb, Tb, La, Ta, U, Para);
    lgnorm = lgnorm + log(Ns);
    
    maxeig = BilayerTraceLTRG(Ta, La, Tb, Lb);
    LnZ(it) = lgnorm + log(maxeig)/2;

%     LnZ(it) = lgnorm;

    fprintf('LnZ: %6f \n', LnZ(it));
    f = @(k) log(1+exp(Para.tau * 2 * it * cos(k * pi)));
    exact(it) = integral(f, 0, 1);
end

% Check result
loglog(beta, abs(LnZ - exact)./exact, '-o', 'linewidth', 2);
xlabel('$\beta$', 'Interpreter', 'latex')
ylabel('$|\delta f/f|$', 'Interpreter', 'latex')

set(gca, 'fontsize', 20, 'fontname', ...
    'times new roman', 'linewidth', 1.5, 'XColor', 'k', 'YColor', 'k')
function [ Ta, Lb, Tb, Ns ] = RG( La, Ta, Lb, Tb, Uab, Para)
% function [ La, Ta, Lb, Tb ] = EvoRG( La, Ta, Lb, Tb, Uab, Para )
D = length(La(:,1));

T = contract(La, 2, Ta, 1);
T = contract(T, 2, Lb, 1, [1,4,2,3]);
T = contract(T, 2, Tb, 1);
T = contract(T, 4, La, 1, [1,2,3,6,4,5]);
T = contract(Uab, [3,4], T, [2,5], [3,1,4,5,2,6]);
T = reshape(T, [D * Para.d^2, D * Para.d^2]);

% //truncate the enlarged bond space by Para.D
[U, S, V, tr, ~] = svdT(T, 'Nkeep', Para.D, 'epsilon', 1e-8);
fprintf('truncation err: %.5e| ', tr);
Dp = length(S);
Ns = norm(diag(S));
Lb = diag(S) / Ns;
U = reshape(U, [D,Para.d,Para.d,Dp]);
U = permute(U, [1,4,2,3]);
V = reshape(V, [Dp,D,Para.d,Para.d]);
Ta = contract(La^(-1), 2, U, 1);
Tb = contract(V, 2, La^(-1), 1, [1,4,2,3]);

end