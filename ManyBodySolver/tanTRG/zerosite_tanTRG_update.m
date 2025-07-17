function [ rho, EE ] = zerosite_tanTRG_update(EnV, rho, T, pos, dir, tau)
fun_prod = @(T) zerosite_prod(T, EnV, pos, dir);
fun_inner = @(T1, T2) tensor_inner_prod(T1, T2);
fun_norm = @(T) tensor_norm(T);
[R, Info] = LanczosExp(fun_prod, T, -tau, 'innerprod', fun_inner, 'norm', fun_norm);
[~, S, ~] = svd(R);
S = diag(S);
S = S/norm(S);
EE = -sum(S.^2 .* log(S.^2));

% updata rho
switch dir
    case '->-'
        rho.A{pos+1} = contract(R, 2, rho.A{pos+1}, 1);
    case '-<-'
        if pos == 2
            rho.A{pos-1} = contract(R, 1, rho.A{pos-1}, 1);
        else
            rho.A{pos-1} = contract(R, 1, rho.A{pos-1}, 2, [2,1,3,4]);
        end
end

end

function [ T ] = zerosite_prod(T, EnV, i, dir)
switch dir
    case '->-'
        T = contract(EnV.A{i}, 3, T, 1);
        T = contract(T, [2,3], EnV.A{i+1}, [2,3]);
    case '-<-'
        T = contract(EnV.A{i-1}, 3, T, 1);
        T = contract(T, [2,3], EnV.A{i}, [2,3]);
end
end

