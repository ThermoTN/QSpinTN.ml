load MBRslt/Rslt.mat
color_list = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]}

load Cm_0T.mat
p(1) = plot(data(:,1), data(:,2), 's', 'color', color_list{1}); hold on
p(2) = plot(Rslt_0T.T, Rslt_0T.Cm * 8.314, '-o', 'linewidth', 2, 'color', color_list{1}); hold on

load Cm_0.05T.mat
plot(data(:,1), data(:,2), 's', 'color', color_list{2}); hold on
p(3) = plot(Rslt_0p05T.T, Rslt_0p05T.Cm * 8.314, '-o', 'linewidth', 2, 'color', color_list{2}); hold on

load Cm_0.1T.mat
plot(data(:,1), data(:,2), 's', 'color', color_list{3}); hold on
p(4) = plot(Rslt_0p1T.T, Rslt_0p1T.Cm * 8.314, '-o', 'linewidth', 2, 'color', color_list{3}); hold on


set(gca, 'fontsize', 20, 'fontname', 'times new roman')
xlabel('$T$ (K)', 'Interpreter', 'latex')
ylabel('$C_m$ (J/mol/K)', 'Interpreter', 'latex')

leg = legend(p, {'Exp.', 'Sim. (0 T)', 'Sim. (0.05 T)', 'Sim. (0.1 T)'}, 'box', 'off');
title(leg, {['$J_{xy}=', num2str(ParaJ(1)*1000, 2), '$ mK, $J_z=', num2str(ParaJ(2)*1000, 2), '$ mK'], ...
    ['$J_{\rm PD}=', num2str(ParaJ(3)*1000, 2), '$ mK, $J_\Gamma=', num2str(ParaJ(4)*1000, 2), '$ mK'], ...
    ['$g_z=', num2str(ParaJ(5),3), '$']}, 'Interpreter', 'latex')
title({'BYB fitting YC3$\times$4'}, 'Interpreter', 'latex')
xlim([0, 0.2])