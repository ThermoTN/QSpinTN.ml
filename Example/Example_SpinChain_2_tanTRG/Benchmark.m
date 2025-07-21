clear all; clc
load Rslt/'Model=XXZtest-Jxy=1-Jz=1-h=[0 0 0.1]-L=12-BC=OBC-D=200.mat'

[En] = GetBenchmarkEDRslt(Para, Rslt.beta);

figure(1)
semilogx(1./Rslt.beta, Rslt.En, '-o'); hold on
plot(1./Rslt.beta, En, '-s')

xlabel('$T$', 'Interpreter', 'latex')
ylabel('$E$', 'Interpreter', 'latex')
legend({['tanTRG, $D=', num2str(Para.MCrit), '$'], 'ED'}, 'Interpreter', 'latex', 'box', 'off', 'location', 'southeast')
set(gca, 'XColor', 'k', 'YColor', 'k', 'fontsize', 20, 'fontname', 'times new roman', 'linewidth', 1.5)

figure(2)
loglog(1./Rslt.beta(2:end), abs(Rslt.En(2:end)-En(2:end)), '-o'); hold on

xlabel('$T$', 'Interpreter', 'latex')
ylabel('$|\delta E|$', 'Interpreter', 'latex')
set(gca, 'XColor', 'k', 'YColor', 'k', 'fontsize', 20, 'fontname', 'times new roman', 'linewidth', 1.5)

function [En] = GetBenchmarkEDRslt(Para, beta_list)
Op = GetLocalSpace(Para.d);
H = AutomataInit(Para);
[H, Ope] = AutomataInit1Site(H, Op, Para);
H = MPO2H(H.A) * exp(H.lgnorm);
[~, D] = eig(H);
D = diag(D);
minD = min(D);
D = D - min(D);
En = beta_list;
for it = 1:length(beta_list)
    En(it) = sum(exp(-beta_list(it) * D) .* D)/sum(exp(-beta_list(it) * D)) + minD;
end
En = En/Para.L;
end