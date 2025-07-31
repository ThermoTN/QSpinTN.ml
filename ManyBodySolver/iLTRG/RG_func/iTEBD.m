function [ LnZ, M ] = iTEBD( Para, T )
% function [ LnZ, M ] = iTEBD( Para )
% Projection: Project the transfer tensors onto the base tensors
%   The base tensors are Ta Tb, with singulr values stored Lama, Lamb.
%   Projectors are Pja, Pjb, return norm factors are the largest singular
%   values in Lama, Lamb (named norma, normb respectively). 
%   projection_num = steps of trotter slices
%   ln(Z)/L = [ln(norma)+ln(normb)]/2.0
%
%      h_ba    |     h_ab     |     h_ba     |      h_ab
%   ----<>----[ ]-----<>-----[ ]-----<>-----[ ]------<>---
%              |              |              |       
%      lama    Ta    lamb     Tb    lama     Ta     lamb
% QMagen 20-March-2021

len = length(datestr(now,'YYYYmmDD_HHMMSS'));

if len == length(Para.TStr_log) && strcmp(Para.WorkingMode, 'ThDQ')
    fileID = 1;
else
    fileID = fopen(['../Tmp/tmp_', Para.TStr, '/', Para.TStr_log, '.log'], 'a');
end

LnZ = zeros(Para.N_max, 1);
M = zeros(Para.N_max, 1);
lgnorm = 0;

[Ta, Tb, La, Lb] = InitId(Para);
[ t_ab, t_ba ] = eval([Para.Trotter_Name, '(Para)']);
[Hab, Hba] = InitHam_trotter(t_ab, t_ba, Para);
[Uab, Uba] = EvoGate(Hab, Hba, Para);
beta = 1./T;
for It = 1:1:Para.N_max
    
    fprintf(fileID, 'It: %6d| beta: %6f| ', It, beta(It));
    
    [Ta, Lb, Tb, Ns] = EvoRG(La, Ta, Lb, Tb, Uab, Para, fileID);
    lgnorm = lgnorm + log(Ns);
    [Tb, La, Ta, Ns] = EvoRG(Lb, Tb, La, Ta, Uba, Para, fileID);
    lgnorm = lgnorm + log(Ns);
    maxeig = BilayerTraceLTRG(Ta, La, Tb, Lb);
    LnZ(It) = lgnorm + log(maxeig)/2;
    
    if norm(Para.Field.h) ~= 0
        M(It) = GetiLTRGM(Ta, La, Tb, Lb, Para);
    else
        M(It) = 0;
    end
    fprintf(fileID, 'LnZ: %6f \n', LnZ(It));
end
 

end

