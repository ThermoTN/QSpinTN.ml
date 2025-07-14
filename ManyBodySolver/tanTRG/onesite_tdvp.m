function [rho, Rslt] = onesite_tdvp(rho, H, Para, Rslt, Op, EnV)

fprintf(Para.fileID, '==== one site TSRG =======================================\n');
while Rslt.beta(end) < Para.beta_max
    Para.tau = GetNowtau(Para, Rslt.beta(end));
    delta = Para.tau * exp(H.lgnorm);
    tic
    for i = 1:1:(Para.L)
        EnV.A{i} = H.A{i};
        [rho, EnV, R] = onesite_tdvp_update(EnV, rho, i, '->-', delta);
        if i ~= Para.L
            [rho, EE(i)] = zerosite_tdvp_update(EnV, rho, R, i, '->-', -delta);
        end
    end
    
    for i = Para.L:-1:1
        EnV.A{i} = H.A{i};
        [rho, EnV, R] = onesite_tdvp_update(EnV, rho, i, '-<-', delta);
        if i ~= 1
            [rho, EE(i-1)] = zerosite_tdvp_update(EnV, rho, R, i, '-<-', -delta);
        end
    end
    
    Rslt.LnZ(end + 1) = 2 * rho.lgnorm / Para.L;
    Rslt.beta(end + 1) = Rslt.beta(end) + 4 * Para.tau;
    [ Rslt ] = caltdvpThDQ( Para, rho, Op, Rslt);
    Rslt.En(end + 1) = GetEn(EnV.A{2}, H, rho) * exp(H.lgnorm) / Para.L;
    Rslt.TrunErr(end + 1) = 0;
    Rslt.EE(end + 1) = max(EE);
    % Obs calculation
    
    % [ Obs ] = GetTSCorr_RealSpace(rho , 1, 5, Op.SOp.Sx, Op.SOp.Sx);
    Rslt = GettanTRGCorr(rho, Rslt, Para, Op);
    savetanTRGRslt(Para, Rslt);
    time = toc();
    fprintf(Para.fileID, 'beta: %f, LnZ: %.16f, En = %.16f\n    Ent. entropy = %f, time = %f s\n', ...
        Rslt.beta(end), Rslt.LnZ(end), Rslt.En(end), max(EE), time);
end

end


