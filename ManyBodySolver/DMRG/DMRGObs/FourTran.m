% function [Sk] = FourTran(Spin, Lx, Ly)
% 2022/4/15 by ezlv
% Spin.xx is the corralate of <SxSx>
close all
Lx = 15;
Ly = 6;
L = Lx * Ly;
% R is a group to save position
% the 1st coml is x position and 2ed coml is y position
R = zeros(L, 2);
a1 = [1, 0];
a2 = [-1/2, -sqrt(3)/2];

for i = 1:L
    cy = mod(i, Ly);
    cx = fix(i/Ly) + 1;
    R(i, 1) = (cx-1)*a1(1) + (cy-1)*a2(1);
    R(i, 2) = (cx-1)*a1(2) + (cy-1)*a2(2);
end

% calculate Sk
Df = 50;
KX = linspace(-2*pi, 2*pi, Df);
KY = linspace(-2*pi, 2*pi, Df);
Sk = zeros(Df, Df);
for kx = 1:Df
    for ky = 1:Df
        for i = 1:L
            for j = 1:L
                if i ~= j
                    Sk(ky, kx) = Sk(ky, kx) + exp(1i * (KX(kx) * (R(i,1) - R(j,1)) + KY(ky) * (R(i,2) - R(j,2))))...
                    * (Spin.xx(min(i,j),max(i,j)) + Spin.yy(min(i,j),max(i,j)) + Spin.zz(min(i,j),max(i,j)))/L;
                else
                    Sk(ky, kx) = Sk(ky, kx) + exp(1i * (KX(kx) * (R(i,1) - R(j,1)) + KY(ky) * (R(i,2) - R(j,2))))...
                    * (3/4)/L;
                end
            end
        end
    end
end
KXf = zeros(Df, Df);
KYf = zeros(Df, Df);
for i = 1:Df
    KXf(i,:) = KX;
    KYf(:,i) = KY;
end
FT = surf(KXf/pi, KYf/pi, real(Sk));
colormap jet
shading interp
colorbar
set(FT, 'EdgeColor', 'none')

% plot 1st BZ
hold on
plot3([-2/3, 2/3], [2/sqrt(3), 2/sqrt(3)], [5, 5], '-', 'color', 'white');
hold on
plot3([2/3, 4/3], [2/sqrt(3), 0], [5, 5], '-', 'color', 'white');
hold on
plot3([2/3, 4/3], [-2/sqrt(3), 0], [5, 5], '-', 'color', 'white');
hold on
plot3([2/3, -2/3], [-2/sqrt(3), -2/sqrt(3)], [5, 5], '-', 'color', 'white');
hold on
plot3([-4/3, -2/3], [0, -2/sqrt(3)], [5, 5], '-', 'color', 'white');
hold on
plot3([-4/3, -2/3], [0, 2/sqrt(3)], [5, 5], '-', 'color', 'white');
xlabel('k_x/\pi', 'FontSize', 10);
ylabel('k_y/\pi', 'FontSize', 10);
% end