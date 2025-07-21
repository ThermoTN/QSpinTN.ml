function [ RG_Rslt ] = GetXTRGRslt( Para )

RG_Rslt.LnZ = [ ];
RG_Rslt.M = [ ];
RG_Rslt.En = [ ];
RG_Rslt.beta = [ ];
RG_Rslt.EE = { };
RG_Rslt.TrunErr = { };
[H, Id, Op] = InitHam(Para);
Op.H = H;
switch Para.ThDQ
    case 'Cm'
        ST_l = [0,1,2,3];
    case 'Chi'
        ST_l = [0,2];
    otherwise
        ST_l = [0];
end


for i = 1:1:length(ST_l)
    ST = ST_l(i);
    Rslt{i} = XTRG_Main( Para, H, Id, Op, ST );
end

for i = 1:1:length(ST_l)
    RG_Rslt.beta = [RG_Rslt.beta; Rslt{i}.beta];
    RG_Rslt.LnZ = [RG_Rslt.LnZ; Rslt{i}.LnZ];
    RG_Rslt.M = [RG_Rslt.M; Rslt{i}.M];
    RG_Rslt.En = [RG_Rslt.En; Rslt{i}.En];
    RG_Rslt.EE = [RG_Rslt.EE, Rslt{i}.EE];
    RG_Rslt.TrunErr = [RG_Rslt.TrunErr, Rslt{i}.TrunErr];
end
RG_Rslt = sort_beta(RG_Rslt);
end

function [ Rslt ] = sort_beta( Rslt )
[Rslt.beta, ord] = sort(Rslt.beta);
Rslt.LnZ = Rslt.LnZ(ord);
Rslt.M = Rslt.M(ord);
Rslt.En = Rslt.En(ord);
Rslt.EE = Rslt.EE(ord);
Rslt.TrunErr = Rslt.TrunErr(ord);
end
