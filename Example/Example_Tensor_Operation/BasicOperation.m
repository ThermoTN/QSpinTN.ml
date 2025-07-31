addpath(genpath('../../ManyBodySolver'))
clear all
clc
%% tensor
T1 = rand([50, 1]);
T2 = rand([50, 50]);
T3 = rand([4, 2, 2]);
T4 = rand([10, 10, 2, 2]);

%% contract
A = rand([10, 10, 2, 2]);
B = rand([10, 10, 2, 2]);
C = contract(A, 2, B, 1);
size(C)

%% permute
C1 = permute(C, [1,4,2,5,3,6]);
size(C1)

%% reshape
sC = size(C1);
D = reshape(C1, [sC(1), sC(2), sC(3) * sC(4), sC(5) * sC(6)]);
size(D)

%% svd
sC = size(C);
T = reshape(C, [prod(sC(1:3)), prod(sC(4:6))]);
[u,S,v] = svd(T); % Note: u * S * v' = T;
sS = size(S);
U = reshape(u, [sC(1), sC(2), sC(3), sS(1)]);
U = permute(U, [1,4,2,3]);
Vd = reshape(v', [sS(2), sC(4), sC(5), sC(6)]);
T = contract(U, 2, S, 1);
T = contract(T, 4, Vd, 1);
norm(reshape(C - T, [numel(C), 1]))/norm(reshape(C, [numel(C), 1]))

%% svdT
sC = size(C);
T = reshape(C, [prod(sC(1:3)), prod(sC(4:6))]);
[u,s,v] = svdT(T); % Note: u * diag(s) * v = T;
S = diag(s);
sS = size(S);
U = reshape(u, [sC(1), sC(2), sC(3), sS(1)]);
U = permute(U, [1,4,2,3]);
Vd = reshape(v, [sS(2), sC(4), sC(5), sC(6)]);
T = contract(U, 2, S, 1);
T = contract(T, 4, Vd, 1);
norm(reshape(C - T, [numel(C), 1]))/norm(reshape(C, [numel(C), 1]))