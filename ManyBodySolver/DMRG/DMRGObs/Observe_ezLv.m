function [Spin] = Observe_ezLv(GS)
% 2022/4/14 by ezlv

d = 2;

[Sx, Sy, Sz, Id] = SpinOp(d);
Sp = Sx + 1i * Sy;
Sm = Sp';
% % calculate correlate <SaSb> of i and j 
% % get a 2D matrix
Spin.xy = Correlate(Sp, Sm, GS);
Spin.zz = Correlate(Sz, Sz, GS);
% Spin.Q = Correlate(Sp^2, Sm^2, GS);
% 
% % calculate the observation of i
% % get a 1D matrix
Spin.x = ObserveSpin(Sx, GS);
% Spin.y = ObserveSpin(Sy, GS);
Spin.z = ObserveSpin(Sz, GS);
% Spin.NM = ObserveSpin(Sp^2, GS);
% % 
% Spin.xydiag = ObserveSpin(Sp * Sm, GS);
% Spin.zzdiag = ObserveSpin(Sz^2, GS);
% Spin.Qdiag = ObserveSpin(Sp^2 * Sm^2, GS);
end

