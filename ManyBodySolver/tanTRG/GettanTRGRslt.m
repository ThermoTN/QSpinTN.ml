function [ Rslt ] = GettanTRGRslt( Para )
[H, Id, Op, ht] = InitHam(Para);
% try
%     hnm = abs(Para.Field.h * ht);
%     spara = rmfield(Para.Model, 'ES');
%     if isfield(spara, 'g')
%         spara = rmfield(spara, 'g');
%     end
%     spara = struct2array(spara);
%     ESp = max(abs([spara, hnm]));
%     H.lgnorm = H.lgnorm - log(ESp);
%     Para.beta_max = Para.beta_max * ESp;
%     Para.beta_c(end) = Para.beta_max;
%     Judge = 1;
%     error('1');
% catch
%     Judge = 0;
% end
% keyboard;
[rho, Op] = SETTN(Para, H, Id, Op, Para.fileID);

[ EnV ] = InitTSRGEnV(rho, H);

Rslt.beta = 2 * Para.tau0;
Rslt.LnZ = 2 * rho.lgnorm/Para.L;
Rslt.M = [];
[ Rslt ] = caltdvpThDQ( Para, rho, Op, Rslt );

Rslt.En = GetEn(EnV.A{2}, H, rho);
Rslt.TrunErr = 0;
Rslt.EE = 0;
% Para.beta_max = 20;
[rho, Rslt, EnV] = twosite_tdvp(rho, H, Para, Rslt, Op, EnV);
[rho, Rslt] = onesite_tdvp(rho, H, Para, Rslt, Op, EnV);
Rslt.EE = Rslt.EE(2:end-1);
Rslt.TrunErr = Rslt.TrunErr(2:end-1);
% if Judge == 1
%     Rslt.beta = Rslt.beta / ESp;
%     Rslt.En = Rslt.En * ESp;
% end
Rslt.rho = rho;
end

