function [rho, Rslt] = onesite_tdvp(rho, H, Para, Rslt, Op, EnV)

fprintf(Para.fileID, '==== one site tdvp =======================================\n');
while Rslt.time_list(end) < Para.tmax
    delta = -1i * Para.tau * exp(H.lgnorm);
    tic
    for i = 1:1:(Para.L)
        EnV.A{i} = H.A{i};
        [rho, EnV, R] = onesite_tdvp_update(EnV, rho, i, '->-', delta, '1t0');
        if i ~= Para.L
            [rho, EE(i)] = zerosite_tdvp_update(EnV, rho, R, i, '->-', -delta);
        end
    end
    
    for i = Para.L:-1:1
        EnV.A{i} = H.A{i};
        [rho, EnV, R] = onesite_tdvp_update(EnV, rho, i, '-<-', delta, '1t0');
        if i ~= 1
            [rho, EE(i-1)] = zerosite_tdvp_update(EnV, rho, R, i, '-<-', -delta);
        end
    end
    
    Rslt.time_list(end + 1) =  Rslt.time_list(end) + 2 * Para.tau;
    [ Obs ] = GetKSpaceObs(rho, Op.rho, Op.QuOp);
    Rslt.Obs(end + 1) = Obs * exp(rho.lgnorm) * exp(1i * Para.E0 * Rslt.time_list(end));
    Rslt.TrunErr(end + 1) = 0;
    Rslt.EE(end + 1) = max(EE);
    time = toc();
    try
        savetdvpRslt(Para, Rslt)
    catch
        if length(Rslt.time_list) == 2
            warning('Can not save result!')
        end
    end
    fprintf(Para.fileID , 't: %.16f\n    Ent. entropy = %f, time = %f s\n', ...
        Rslt.time_list(end), max(EE), time);
    
    
end


