% 2022/4/15 by ezlv
% Plot Spin x and Spin y in 2D Map
close all
load('Spin_4T.mat');
Lx = 15;
Ly = 6;
L = Lx * Ly;
[LocX, LocY] = meshgrid(1:Lx, 1:Ly);
SpinX = zeros(Ly, Lx);
SpinY = zeros(Ly, Lx);
SpinZ = zeros(Ly, Lx);
a1 = [1, 0];
a2 = [-1/2, -sqrt(3)/2];

for i = 1:Lx
    for j = 1:Ly
        LocX(j, i) = (i-1)*a1(1) + (j-1)*a2(1);
        LocY(j, i) = (i-1)*a1(2) + (j-1)*a2(2);
        SpinX(j, i) = real(Spin.x((i-1)*Ly + j));
        SpinY(j, i) = real(-Spin.y((i-1)*Ly + j));
        SpinZ(j, i) = real(Spin.z((i-1)*Ly + j));
    end
end
% plot meshgrid
for i = 1:Lx
    plot([LocX(1,i), LocX(end,i)], [LocY(1,i), LocY(end,i)], 'k-');
    hold on
end
for i = 1:Ly
    plot([LocX(i,1), LocX(i,end)], [LocY(i,1), LocY(i,end)], 'k-');
    hold on
end
for i = 1:Ly
    plot([LocX(Ly,i), LocX(Ly-i+1,1)], [LocY(Ly,i), LocY(Ly-i+1,1)], 'k-');
    hold on    
end
for i = Ly+1:Lx
    plot([LocX(Ly,i), LocX(1,i-Ly+1)], [LocY(Ly,i), LocY(1,i-Ly+1)], 'k-')
    hold on
end
for i = Lx+1:Lx+Ly-1
    plot([LocX(1,i-Ly+1), LocX(Lx+Ly-i,Lx)], [LocY(1,i-Ly+1), LocY(Lx+Ly-i,Lx)], 'k-')
    hold on
end

for i = 1:Lx
    for j = 1:Ly
        LocX1(j, i) = LocX(j, i) - SpinX(j, i)/2;
        LocY1(j, i) = LocY(j, i) - SpinY(j, i)/2;
    end
end

hold on
q = quiver(LocX1, LocY1, SpinX, SpinY, 0.4, 'red', 'LineWidth', 1.5);
axis equal
axis([-(Ly-1)/2, Lx - 1, fix(-Ly * sqrt(3)/2) , 1]);

hold on
LocX = reshape(LocX, [1, Lx*Ly]);
LocY = reshape(LocY, [1, Lx*Ly]);
SpinZ = reshape(SpinZ, [1, Lx*Ly]);

X = [];
Y = [];
Z = [];
for i = 1:L
    if SpinZ(i) > 0
    X(end+1) = LocX(i);
    Y(end+1) = LocY(i);
    Z(end+1) = 2*SpinZ(i).^2*2000;
    end
end
s1 = scatter(X, Y, Z, 'blue', 'filled');
hold on
%s = scatter(X, Y, 2000, 'yellow');
hold on
X = [];
Y = [];
Z = [];
for i = 1:L
    if SpinZ(i) < 0
    X(end+1) = LocX(i);
    Y(end+1) = LocY(i);
    Z(end+1) = 2*SpinZ(i).^2*2000;
    end
end
s2 = scatter(X, Y, Z, 'yellow', 'filled');
legend([q, s1, s2], 'Sxy', 'Sz>0', 'Sz<0', 'Box', 'off', 'NumColumns',3);
xlabel('L_x');
ylabel('L_y');
title('YbMgGaO_4 GS');