function [ rho, Op ] = SETTN( Para, H, Id, Op, fileID )
% function [ rho, Op ] = SETTN( Para, H, Id, Op )
% Series-Expansion Thermal Tensor Network algorithm to get 
% the density matrix rho at inverse temperature Para.tau.
% rho.A is the tensor that make up the MPO of density martrix.
% rho.lgnorm is the log norm of density martrix.
% Yuan Gao@buaa 2020.12.07
% mail: 17231064@buaa.edu.cn

fprintf(fileID, '===================================================================\n');
fprintf(fileID, 'SETTN:\n');


tau = Para.tau0;
step_max = Para.SETTN_init_step_max;

rho = Id;
Hn = H;
switch Para.ManyBodySolver
    case 'tanTRG'
        Para.MCrit = floor(Para.MCrit/4);
end
for i = 1:1:step_max
    % //rho = rho + (-tau)^i/i! H^i
    sig = sign((-1)^i);
    fec1 = rho.lgnorm;
    fec2 = Hn.lgnorm + i * log(tau) - sum(log(1:1:i));
    
    if fec1 > fec2
        [rho.A, nm, TE, EE] = VariSumMPO(Para, rho.A, Hn.A, [1, sig * exp(fec2-fec1)]);
        rho.lgnorm = fec1 + log(nm);
    else
        [rho.A, nm, TE, EE] = VariSumMPO(Para, rho.A, Hn.A, [exp(fec1-fec2), sig]);
        rho.lgnorm = fec2 + log(nm);
    end
    
    fprintf(fileID, '    %d:        Truncation err of VarSum: %g\n', i, max(TE));
%     fprintf('%d, %.5f, %.5f\n', i, fec1/log(10), fec2/log(10));
    if isfield(Para, 'TSRGStol')
        tol = log(Para.TSRGStol)/log(10);
    else
        tol = -14;
    end
    
    if abs(rho.lgnorm-fec1) < 10^(tol)
        break;
    end
    % //Hn = H * Hn (H^(i+1))
    [Hn.A, nm, TE, EE] = VariProdMPO(Para, Hn.A, H.A);
    Hn.lgnorm = H.lgnorm + Hn.lgnorm + log(nm);
end
fprintf(fileID, '-------------------------------------------------------------------\n');
end
