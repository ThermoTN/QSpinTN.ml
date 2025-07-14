function [ rho, Rslt ] = ExpEvo_conj( Para, rho, Op, fileID )
% function[ rho ] = ExpEvo_conj( Para, rho, H, Op )
% Exponential Thermal Tensor Network algorithm to calculate the partition
% function.
% Yuan Gao@buaa 2020.12.07
% mail: 17231064@buaa.edu.cn
fprintf(fileID, '===================================================================\n');
fprintf(fileID, 'XTRG:\n');

It_max = Para.It;

LnZ = zeros(It_max, 1);
Cm = zeros(It_max, 1);
M = zeros(It_max, 1);
En = zeros(It_max, 1);
beta = zeros(It_max, 1);

global EVOFLAG;

for It = 1:1:It_max
    if EVOFLAG == 1
        tic
    end
    Para.MCrit = Para.D_list(It);
    crho = rho;
    crho.A{1} = conj(rho.A{1});
    crho.A{1} = permute(crho.A{1}, [1,3,2]);
    crho.A{end} = conj(rho.A{end});
    crho.A{end} = permute(crho.A{end}, [1,3,2]);
    for i = 2:1:length(rho.A)-1
        crho.A{i} = conj(crho.A{i});
        crho.A{i} = permute(crho.A{i}, [1,2,4,3]);
    end
    
    % change the max interation step for rho-rho process
    Para.VariProd_step_max = Para.sweep_step;
    [rho.A, nm, TE, ee] = VariProdMPO(Para, rho.A, crho.A);
    rho.lgnorm = (2 * rho.lgnorm + log(nm));
    if EVOFLAG == 1
        toc
    end
    
    
    
    LnZ(It) = 2 * rho.lgnorm / Para.L;
    beta(It) = Para.tau0 * 2 ^ (It+1);
    ThDQ = calThDQXTRG(Para, rho, Op);
    M(It) = ThDQ.M / Para.L;
    En(It) = ThDQ.En / Para.L;
    EE{It} = ee;
    TrunErr{It} = TE;
    TE = max(TE);
    fprintf(fileID, '    %d:    beta=%3.6f, LnZ=%.16f, En=%2.15f\n', It, beta(It), LnZ(It), En(It));
    fprintf(fileID, '        Truncation err %g\n        Entanglement Entropy %.6f\n', TE, ee(floor(Para.L/2)+1));
    
end
fprintf(fileID, '===================================================================\n');
Rslt.LnZ = LnZ;
Rslt.Cm = Cm;
Rslt.M = M;
Rslt.En = En;
Rslt.beta = beta;
Rslt.EE = EE;
Rslt.TrunErr = TrunErr;
end

