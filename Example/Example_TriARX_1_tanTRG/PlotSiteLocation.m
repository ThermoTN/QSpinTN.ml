load Rslt/'Model=TLARX-J1xy=1-J1z=1-JPD=0.5-JGamma=0-J2xy=0-J2z=0-h=[0 0 0.1]-L=12-Lx=4-Ly=3-BCX=OBC-BCY=PBC-D=200.mat'
vec = GetPosi(Para, 1:12);

scatter(vec(:,1), vec(:,2), 50, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on

xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
set(gca, 'XColor', 'k', 'YColor', 'k', 'fontsize', 20, 'fontname', 'times new roman', 'linewidth', 1.5)
axis equal

for it = 1:12
    scatter(vec(it,1), vec(it,2), 50, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); hold on
    pause(0.2)
end

Para.Geo.BCY = 'OBC';
[ Bond_info, ~, ~ ] = Bond_Info_TLARX( Para );
for it = 1:length(Bond_info(:,1))
    vx = [vec(Bond_info{it, 1}, 1), vec(Bond_info{it, 2}, 1)];
    vy = [vec(Bond_info{it, 1}, 2), vec(Bond_info{it, 2}, 2)];
    switch Bond_info{it, 3}
        case 'nn-1'
            plot(vx, vy, '--', 'color', 'r', 'linewidth', 2)
        case 'nn-2'
            plot(vx, vy, '--', 'color', 'g', 'linewidth', 2)
        case 'nn-3'
            plot(vx, vy, '--', 'color', 'b', 'linewidth', 2)
    end
end