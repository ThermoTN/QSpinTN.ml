function [ rho, Rslt, EnV ] = twosite_tdvp(rho, H, Para, Rslt, Op, EnV)

fprintf(Para.fileID, '==== two site tdvp =======================================\n');
while Rslt.time_list(end) < Para.tmax
    delta = -1i * Para.tau * exp(H.lgnorm);
    tic
    for i = 1:1:(Para.L-1)
        EnV.A{i+1} = H.A{i+1};
        [rho, EnV, TE(i), EE(i)] = twosite_tdvp_update(EnV, rho, i, '->-', delta, Para);
        if i ~= Para.L-1
            [rho, EnV, ~] = onesite_tdvp_update(EnV, rho, i+1, '->-', -delta, '2t1');
        end
    end
    
    for i = Para.L:-1:2
        EnV.A{i-1} = H.A{i-1};
        [rho, EnV, TE(i-1), EE(i-1)] = twosite_tdvp_update(EnV, rho, i, '-<-', delta, Para);
        if i ~= 2
            [rho, EnV, ~] = onesite_tdvp_update(EnV, rho, i-1, '-<-', -delta, '2t1');
        end
    end
    
    Rslt.time_list(end + 1) =  Rslt.time_list(end) + 2 * Para.tau;
    [ Obs ] = GetKSpaceObs(rho, Op.rho, Op.QuOp);
    Rslt.Obs(end + 1) = Obs * exp(rho.lgnorm) * exp(1i * Para.E0 * Rslt.time_list(end));
    Rslt.TrunErr(end + 1) = max(TE);
    Rslt.EE(end + 1) = max(EE);
    time = toc();
    D = getBondDim(rho);
    fprintf(Para.fileID , 't: %.16f\n    Ent. entropy = %f, Truncation error = %g, D = %d, time = %f s\n', ...
        Rslt.time_list(end), max(EE), max(TE), D, time);
    
    try
        savetdvpRslt(Para, Rslt)
    catch
        if length(Rslt.time_list) == 2
            warning('Can not save result!')
        end
    end
    
    if Rslt.time_list(end) > Para.t_switch && ~Para.Full2
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
