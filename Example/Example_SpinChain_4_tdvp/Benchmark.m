clear all; clc
load Rslt_DY/'Model=XXZtest-Jxy=1-Jz=1-h=[0 0 0.1]-L=10-BC=OBC-D=100-k=3.14-Op=Sz.mat'

[Obs] = GetBenchmarkEDRslt(Para, Rslt.time_list);

figure(1)
plot(1./Rslt.time_list, real(Rslt.Obs), '-o'); hold on
plot(1./Rslt.time_list, real(Obs), '-s')

xlabel('$t$', 'Interpreter', 'latex')
ylabel('${\rm Re}[g^{zz}(\pi, t)]$', 'Interpreter', 'latex')
legend({['tdvp, $D=', num2str(Para.MCrit), '$'], 'ED'}, 'Interpreter', 'latex', 'box', 'off', 'location', 'southeast')
set(gca, 'XColor', 'k', 'YColor', 'k', 'fontsize', 20, 'fontname', 'times new roman', 'linewidth', 1.5)

figure(2)
semilogy(1./Rslt.time_list, abs(Rslt.Obs-Obs), '-o'); hold on

xlabel('$t$', 'Interpreter', 'latex')
ylabel('$|\delta g^{zz}(\pi, t)|$', 'Interpreter', 'latex')
set(gca, 'XColor', 'k', 'YColor', 'k', 'fontsize', 20, 'fontname', 'times new roman', 'linewidth', 1.5)

function [Obs] = GetBenchmarkEDRslt(Para, time_list)
Op = GetLocalSpace(Para.d);
H = AutomataInit(Para);
[H, Ope] = AutomataInit1Site(H, Op, Para);
H = MPO2H(H.A) * exp(H.lgnorm);
[V, D] = eig(H);
GS = V(:,1);

QOp = 0;
Sz = [-0.5, 0; 0, 0.5];
for it = 1:Para.L
    QOp = QOp + exp(1i * Para.DY.k * it) * kron(eye(Para.d^(Para.L-it)), kron(Sz, eye(Para.d^(it-1))));
end
QOp = QOp/sqrt(Para.L);
QOp1 = AutomataInit_MPO(Para.QuIntr, Para, Para.L);
U = expm(-1i * H * time_list(2));
for it = 1:length(time_list)
    if it == 1
        A = eye(Para.d^Para.L);
    else
        A = A * U;
    end
    Obs(it) = GS' * QOp' * A * QOp * GS * exp(1i * time_list(it) * D(1,1));
end
end