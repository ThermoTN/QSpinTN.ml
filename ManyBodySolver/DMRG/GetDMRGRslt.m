function [Rslt, T] = GetDMRGRslt(Para)
% The optional state output also allows measurement when saveMPS is false.

[ H, ~, Ope ] = InitHam(Para);
[T, E0, EE, envDir, info] = DMRGMain(Para, H);
Rslt.DMRGInfo = info;
if ~isempty(envDir), Rslt.EnvDir = envDir; end
if norm(Para.Field.h) ~= 0
    M = Observe( T, Ope.Sm, Para )/Para.L;
else
    M = 0;
end
[Op.Sx, Op.Sy, Op.Sz, Op.Id] = SpinOp(Para.d);

if isfield(Para,'saveMPS') && Para.saveMPS
    Rslt.T = T;
    Rslt.H = H;
end
Rslt.E0 = E0;
Rslt.M = M;
Rslt.Op = Op;
Rslt.EE = EE;
% [Rslt.Obs] = Observe_ezLv(T.A);

end

