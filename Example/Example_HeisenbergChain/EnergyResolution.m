clear all


close all
set(gcf, 'Position', [1 452 1036 420])

axes('Position', [0.08, 0.15, 0.35, 0.7])
t = -1:0.001:1;
PZ = parzenwin(length(t));
plot(t, PZ, 'linewidth', 3)


set(gca, 'fontsize', 22, 'fontname', 'times new roman', ...
    'XColor', 'k', 'YColor', 'k', 'linewidth', 2)
hold off
xlabel('$t$', 'Fontsize', 22, 'Interpreter', 'latex')
ylabel('$W(t)$', 'Fontsize', 22, 'Interpreter', 'latex')


axes('Position', [0.58, 0.15, 0.35, 0.7])
omega = -10:0.001:10;
plot(omega, FP(omega), 'linewidth', 3); hold on
a = max(FP(omega));

plot([-4, 4], [a,a]/2,  '--', 'linewidth', 3)
set(gca, 'fontsize', 22, 'fontname', 'times new roman', ...
    'XColor', 'k', 'YColor', 'k', 'linewidth', 2)
hold off
xlabel('$\omega$', 'Fontsize', 22, 'Interpreter', 'latex')
ylabel('$F[W(t)](\omega)$', 'Fontsize', 22, 'Interpreter', 'latex')
text(-0.1, 0.13, '$\varepsilon$', 'Interpreter', 'latex', 'Fontsize', 22)

axes('Position', [0,0,1,1])
set(gca, 'XColor', 'none', 'YColor', 'none', 'Color', 'none')
sh = 0.08;
text(0.08-0.07, 0.85, '(a)', 'Fontsize', 22, 'Fontname', 'times new roman')
text(0.58-0.08, 0.85, '(b)', 'Fontsize', 22, 'Fontname', 'times new roman')

function data = FP(omega)
data = real(6 * sqrt(2/pi) * exp(-1i * omega) .*(-1 + exp(1i * omega/2)).^4./omega.^4);
end