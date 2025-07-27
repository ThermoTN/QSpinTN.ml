clear all
k_list = 0:0.0625:1;

ll = 1:101;
omega = 0:0.001:4;
for it = 1:length(k_list)
    load(['Rslt_DY/Model=XXZtest-Jxy=1-Jz=1-h=[0 0 0]-L=64-BC=OBC-D=200-k=', ...
        mat2str(k_list(it) * pi, 3), '-Op=Sz.mat']);
    for itw = 1:length(omega)
        [SObs(it, itw), ChiObs(it, itw)] = Getomega(Rslt.Obs(ll), omega(itw), Rslt.time_list(ll), 'on');
    end
end

SObs = [SObs;SObs(end-1:-1:1,:)] * 3;
ChiObs = [ChiObs;ChiObs(end-1:-1:1,:)] * 3;
kp = 0:0.0625:2;
[xms, yms] = meshgrid(kp, omega);
contourf(xms, yms, SObs', 1000, 'linecolor', 'none')
c = colorbar;
colormap('jet')
c.LineWidth = 1.5;
caxis([0, 20]);
xlabel('$k/\pi$', 'Interpreter', 'latex')
ylabel('$\omega$', 'Interpreter', 'latex')
set(gca, 'fontsize', 22, 'fontname', 'times new roman', 'linewidth', 1.5)
title({'$S(k,\omega)$ of 1D Heisenberg Chain', ['$D=200$, $L=64$, $\varepsilon=', num2str(8/Rslt.time_list(ll(end))),'$']}, 'Interpreter', 'latex')


hold on

k = 0:0.001:2;

plot(k, pi * abs(sin(k * pi))/2, 'color', 'w', 'linewidth', 2); 
plot(k, pi * abs(sin(k * pi/2)), 'color', 'w', 'linewidth', 2); 
hold off
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
