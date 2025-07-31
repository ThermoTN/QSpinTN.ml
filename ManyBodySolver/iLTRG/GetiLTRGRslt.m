function [ Rslt ] = GetiLTRGRslt( Para, ThDQ )

[Rslt.beta, Rslt.T] = Tem_cal(Para);  % get temperature points

% //Compute specific heat
if strcmp(ThDQ, 'Cm') || strcmp(ThDQ, 'Cm&Chi')
    [LnZ_l, Ml] = iTEBD(Para, Rslt.T);              % call iLTRG to compute lnZ
    Rslt.Cm = Cm_cal(LnZ_l, Rslt.beta, Para);
    Rslt.Cm = Rslt.Cm(2:1:(end-1));
    Rslt.M = Ml(2:1:(end-1));
else
    Rslt.Cm = 0;
end
Rslt.LnZ = LnZ_l(2:1:(end-1));
Rslt.beta = Rslt.beta(2:1:(end-1));
Rslt.T = Rslt.T(2:1:(end-1));

Tp = Rslt.T;
T_min = min(Rslt.T) * 0.99;
T_max = max(Rslt.T) * 0.99;
Delta = (log(T_max)-log(T_min))/Para.InterNum;
Rslt.T = exp(log(T_max):-Delta:log(T_min));

if strcmp(ThDQ, 'Cm') || strcmp(ThDQ, 'Cm&Chi')
    Rslt.Cm = interp1(Tp, Rslt.Cm, Rslt.T);
    Rslt.Cm = transpose(Rslt.Cm);
end

if norm(Para.Field.h) ~= 0 && (strcmp(ThDQ, 'Chi') || strcmp(ThDQ, 'Cm&Chi'))
    Rslt.M = interp1(Tp, Rslt.M, Rslt.T);
    Rslt.M = transpose(Rslt.M);
else
    Rslt.M = 0;
end
Rslt.T = transpose(Rslt.T);
Rslt.beta = 1 ./ Rslt.T;
end

