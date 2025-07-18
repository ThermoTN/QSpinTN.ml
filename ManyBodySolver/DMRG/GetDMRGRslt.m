function [Rslt] = GetDMRGRslt(Para)

[ H, Id, Ope ] = InitHam(Para);
[T, E0] = DMRGMain(Para, H);
if norm(Para.Field.h) ~= 0
    M = Observe( T, Ope.Sm, Para )/Para.L;
else
    M = 0;
end
[Op.Sx, Op.Sy, Op.Sz, Op.Id] = SpinOp(Para.d);

try
    if Para.saveMPS
        Rslt.T = T;
        Rslt.H = H;
    end
end
Rslt.E0 = E0;
Rslt.M = M;
Rslt.Op = Op;

% [Rslt.Obs] = Observe_ezLv(T.A);

end

