load Rslt/'Model=TLARX-J1xy=1-J1z=1-JPD=0.5-JGamma=0-J2xy=0-J2z=0-h=[0 0 0.1]-L=12-Lx=4-Ly=3-BCX=OBC-BCY=PBC-D=200.mat'

[Cm, beta] = ThDQ_func(Rslt);

semilogx(1./beta, Cm, '-o', 'linewidth', 2)

xlim([0, 10])

xlabel('$T$', 'Interpreter', 'latex')
ylabel('$C_{\rm m}$', 'Interpreter', 'latex')

set(gca, 'XColor', 'k', 'YColor', 'k', 'fontsize', 20, 'fontname', 'times new roman', 'linewidth', 1.5)