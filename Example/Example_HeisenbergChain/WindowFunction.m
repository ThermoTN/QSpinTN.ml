clear all
clc
figure(1)
ll = 1:51;

[omega, SObs] = PlotLine(ll, 'on');
plot(omega, SObs, 'linewidth', 2, 'color', 'b'); hold on
[omega, SObs] = PlotLine(ll, 'off');
plot(omega, SObs, '-.', 'linewidth', 2, 'color', 'b'); hold on

ll = 1:101;
[omega, SObs] = PlotLine(ll, 'on');
plot(omega, SObs, 'linewidth', 2, 'color', 'r'); hold on
[omega, SObs] = PlotLine(ll, 'off');
plot(omega, SObs, '-.', 'linewidth', 2, 'color', 'r'); hold on

xlabel('$k/\pi$', 'Interpreter', 'latex')
ylabel('$\omega$', 'Interpreter', 'latex')
set(gca, 'fontsize', 22, 'fontname', 'times new roman', 'linewidth', 1.5)

hold off
ylim([0, 30])
legend({'$t_{\rm max} = 10$, Parzen', 'Staircase', ...
    '$t_{\rm max} = 20$, Parzen', 'Staircase'}, 'box', 'off', 'Interpreter', 'latex')
title('$S(\pi, \omega)$', 'Interpreter', 'latex');

figure(2)

x = -1:0.001:1;
plot([-5, x, 5], [0, parzenwin(length(x))', 0], '-', 'linewidth', 2); hold on

plot([-5, -1, -1, 1, 1, 5], [0,0,1,1,0,0], '-', 'linewidth', 2); hold on

xlim([-3, 3])
xlabel('$t/t_{\rm max}$', 'Interpreter', 'latex')
ylabel('$W(t/t_{\rm max})$', 'Interpreter', 'latex')
legend({'Parzen', 'Staircase'}, 'box', 'off')
set(gca, 'fontsize', 22, 'fontname', 'times new roman', 'linewidth', 1.5)
hold off

function [omega, SObs] = PlotLine(ll, win)
omega = 0:0.001:4;
load(['Rslt_DY/Model=XXZtest-Jxy=1-Jz=1-h=[0 0 0]-L=64-BC=OBC-D=200-k=', ...
    mat2str(pi, 3), '-Op=Sz.mat']);
for itw = 1:length(omega)
    [SObs(itw), ~] = Getomega(Rslt.Obs(ll), omega(itw), Rslt.time_list(ll), win);
end
end

function [OmegaSObs, OmegaChiObs] = Getomega(ObsData, omega, time_list, win)
ObsData = reshape(ObsData, [1, numel(ObsData)]);
window = parzenwin(2 * length(time_list) - 1);
window = window(length(time_list):1:end);

tran = 2 * exp(1i * omega * time_list);
tranp = 2 * sin(omega * time_list);
switch win
    case 'on'
        tran = tran .* window';
        tranp = tranp .* window';
    case 'off'
        
end
OmegaSObs = trapz(time_list, real(tran .* ObsData));
OmegaChiObs = -trapz(time_list, tranp .* imag(ObsData));
end
