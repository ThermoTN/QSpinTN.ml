function [ Obs ] = GetTSCorr_RealSpace(rho , loc1, loc2, Op1, Op2)
rhod = rho;
switch loc1
    case {1, length(rho.A)}
        rho.A{loc1} = contract(rho.A{loc1}, 2, Op1, 2, [1,3,2]);
    otherwise
        rho.A{loc1} = contract(rho.A{loc1}, 3, Op1, 2, [1,2,4,3]);
end

switch loc2
    case {1, length(rho.A)}
        rho.A{loc2} = contract(rho.A{loc2}, 2, Op2, 2, [1,3,2]);
    otherwise
        rho.A{loc2} = contract(rho.A{loc2}, 3, Op2, 2, [1,2,4,3]);
end
T = contract(rho.A{1}, [2,3], conj(rhod.A{1}), [2,3]);
if max([loc1, loc2]) < length(rho.A)
    for It = 2:max([loc1, loc2])
        T = contract(T, 1, rho.A{It}, 1);
        T = contract(T, [1,3,4], conj(rhod.A{It}), [1,3,4]);
    end
    Obs = trace(T);
else
    for It = 2:length(rho.A) - 1
        T = contract(T, 1, rho.A{It}, 1);
        T = contract(T, [1,3,4], conj(rhod.A{It}), [1,3,4]);
    end
    T = contract(T, 1, rho.A{end}, 1);
    Obs = contract(T, [1,2,3], conj(rhod.A{end}), [1,2,3]);
end
end

