function [ ED_Rslt ] = GetEDRslt( Para )

Op = GetLocalSpace(Para.d);
H = AutomataInit(Para);
[H, Ope] = AutomataInit1Site(H, Op, Para);
H = MPO2H(H.A) * exp(H.lgnorm);
beta_list = Para.beta_list;

ED_Rslt.H = H;
ThDQ = Para.ThDQ;
ED_Rslt.Op = Op;
switch ThDQ
    case 'Cm'
        [T, C, E, F] = ED_Cm(H, beta_list);
        C = C ./ Para.L;
        ED_Rslt.T = T;
        ED_Rslt.beta = 1./T;
        ED_Rslt.Cm = C;
        ED_Rslt.M = 0;
        ED_Rslt.En = E ./ Para.L;
        ED_Rslt.F = F ./ Para.L;
    case 'Chi'
        M = MPO2H(Ope.Sm.A);
        ED_Rslt.Mop = M;
        [T, M, ~] = ED_chi(H, M, norm(Para.Field.h), beta_list);
        M = M ./ Para.L;
        ED_Rslt.T = T;
        ED_Rslt.beta = 1./T;
        ED_Rslt.M = M;
        ED_Rslt.Cm = 0;
    otherwise
        if norm(Para.Field.h) ~= 0
            M = MPO2H(Ope.Sm.A);
            ED_Rslt.Mop = M;
        end
end
end

