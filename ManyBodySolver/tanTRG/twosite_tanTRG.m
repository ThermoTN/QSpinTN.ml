function [ rho, Rslt, EnV ] = twosite_tanTRG(rho, H, Para, Rslt, Op, EnV)

fprintf(Para.fileID, '==== two site tanTRG =======================================\n');
while Rslt.beta(end) < Para.beta_max
    Para.tau = GetNowtau(Para, Rslt.beta(end));
    delta = Para.tau * exp(H.lgnorm);
    tic
    for i = 1:1:(Para.L-1)
        EnV.A{i+1} = H.A{i+1};
        [rho, EnV, TE(i), EE(i)] = twosite_tanTRG_update(EnV, rho, i, '->-', delta, Para);
        if i ~= Para.L-1
            [rho, EnV, ~] = onesite_tanTRG_update(EnV, rho, i+1, '->-', -delta);
        end
    end
    
    for i = Para.L:-1:2
        EnV.A{i-1} = H.A{i-1};
        [rho, EnV, TE(i-1), EE(i-1)] = twosite_tanTRG_update(EnV, rho, i, '-<-', delta, Para);
        if i ~= 2
            [rho, EnV, ~] = onesite_tanTRG_update(EnV, rho, i-1, '-<-', -delta);
        end
    end
    
    Rslt.LnZ(end + 1) = 2 * rho.lgnorm / Para.L;
    Rslt.beta(end + 1) =  Rslt.beta(end) + 4 * Para.tau;
    [ Rslt ] = caltdvpThDQ( Para, rho, Op, Rslt);
    Rslt.En(end + 1) = GetEn(EnV.A{2}, H, rho) * exp(H.lgnorm)/ Para.L;
    Rslt.TrunErr(end + 1) = max(TE);
    Rslt.EE(end + 1) = max(EE);
    time = toc();
    D = getBondDim(rho);
    fprintf(Para.fileID , 'beta: %.16f, LnZ: %.16f, En = %.16f\n    Ent. entropy = %f, Truncation error = %g, D = %d, time = %f s\n', ...
        Rslt.beta(end), Rslt.LnZ(end), Rslt.En(end), max(EE), max(TE), D, time);
    try
        Rslt = GettanTRGCorr(rho, Rslt, Para, Op);
        savetanTRGRslt(Para, Rslt);
    catch
        if length(Rslt.beta) == 2
            warning('Unable to save the results / calculate observable!\n');
        end
    end
    
    if Rslt.beta(end) > Para.beta_switch && ~Para.Full2
        break
    end
end
end

function D = getBondDim(rho)
D = 0;
for It = 1:length(rho.A)
    D = max(D, size(rho.A{It}, 1));
end
end
