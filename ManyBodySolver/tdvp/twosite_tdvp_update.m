function [ rho, EnV, TE, EE ] = twosite_tdvp_update(EnV, rho, pos, dir, tau, Para)
D_max = Para.MCrit;
switch dir
    case '->-'
        if pos == 1
            T = contract(rho.A{pos}, 1, rho.A{pos+1}, 1, [2,1,3]);
        elseif pos == length(EnV.A) - 1
            T = contract(rho.A{pos}, 2, rho.A{pos+1}, 1);
        else
            T = contract(rho.A{pos}, 2, rho.A{pos+1}, 1, [1,3,2,4]);
        end
    case '-<-'
        if pos == length(EnV.A)
            T = contract(rho.A{pos-1}, 2, rho.A{pos}, 1);
        elseif pos == 2
            T = contract(rho.A{pos-1}, 1, rho.A{pos}, 1, [2,1,3]);
        else
            T = contract(rho.A{pos-1}, 2, rho.A{pos}, 1, [1,3,2,4]);
        end
end

fun_prod = @(T) twosite_prod(T, EnV, pos, dir);
fun_inner = @(T1, T2) tensor_inner_prod(T1, T2);
fun_norm = @(T) tensor_norm(T);
[T, Info] = LanczosExp(fun_prod, T, tau, 'innerprod', fun_inner, 'norm', fun_norm);
% [T, Info] = LanczosExp(fun_prod, T, -tau, 'innerprod', fun_inner, 'norm', fun_norm);

len = length(rho.A);
if pos == 1 || (pos == 2 && strcmp(dir, '-<-'))
    T = permute(T, [2,1,3]);
    Tsize = size(T);
    T = reshape(T, [prod(Tsize(1)), prod(Tsize(2:3))]);
elseif (pos == (len-1) && strcmp(dir, '->-')) || pos == len
    Tsize = size(T);
    T = reshape(T, [prod(Tsize(1:2)), prod(Tsize(3))]);
else
    T = permute(T, [1,3,2,4]);
    Tsize = size(T);
    T = reshape(T, [prod(Tsize(1:2)), prod(Tsize(3:4))]);
end

[U, S, V, TE, EE] = svdT(T, 'Nkeep', D_max, 'epsilon', Para.Stol);
Usize = size(U);
Vsize = size(V);
if pos == 1 || (pos == 2 && strcmp(dir, '-<-'))
    U = reshape(U, [Tsize(1), Usize(2)]);
    U = permute(U, [2,1]);
    V = reshape(V, [Vsize(1), Tsize(2), Tsize(3)]);
elseif (pos == (len-1) && strcmp(dir, '->-')) || pos == len
    U = reshape(U, [Tsize(1), Tsize(2), Usize(2)]);
    U = permute(U, [1,3,2]);
    V = reshape(V, [Vsize(1), Tsize(3)]);
else
    U = reshape(U, [Tsize(1), Tsize(2), Usize(2)]);
    U = permute(U, [1,3,2]);
    V = reshape(V, [Vsize(1), Tsize(3), Tsize(4)]);
end
% fprintf('Trunction Error: %g, min S: %g\n', TE, min(S));
Ns = norm(S);
S = S/Ns;
S = diag(S);
switch dir
    case '->-'
        rho.A{pos} = U;
        rho.A{pos+1} = contract(S, 2, V, 1);
    case '-<-'
        if pos-1 == 1
            rho.A{pos-1} = contract(U, 1, S, 1, [2,1]);
        else
            rho.A{pos-1} = contract(U, 2, S, 1, [1,3,2]);
        end
        rho.A{pos} = V;
end
log(Ns);
rho.lgnorm = rho.lgnorm + log(Ns);

% update EnV
switch dir
    case '->-'
        if pos == 1
            EnV.A{pos} = contract(rho.A{pos}, 2, EnV.A{pos}, 3);
            EnV.A{pos} = contract(EnV.A{pos}, 3, conj(rho.A{pos}), 2, [3,2,1]);
        else
            EnV.A{pos} = contract(EnV.A{pos-1}, 2, EnV.A{pos}, 1);
            EnV.A{pos} = contract(EnV.A{pos}, [1,4], conj(rho.A{pos}), [1,3]);
            EnV.A{pos} = contract(EnV.A{pos}, [1,3], rho.A{pos}, [1,3], [2,1,3]);
        end
    case '-<-'
        if pos == length(EnV.A)
            EnV.A{pos} = contract(rho.A{pos}, 2, EnV.A{pos}, 3);
            EnV.A{pos} = contract(EnV.A{pos}, 3, conj(rho.A{pos}), 2, [3,2,1]);
        else
            EnV.A{pos} = contract(EnV.A{pos+1}, 2, EnV.A{pos}, 2);
            EnV.A{pos} = contract(EnV.A{pos}, [1,4], conj(rho.A{pos}), [2,3]);
            EnV.A{pos} = contract(EnV.A{pos}, [1,3], rho.A{pos}, [2,3], [2,1,3]);
        end
end
end

function [ T ] = twosite_prod(T, EnV, i, dir)
switch dir
    case '->-'
        if i == 1
            T = contract(T, 1, EnV.A{i+2}, 3);
            T = contract(T, 1, EnV.A{i}, 3);
            T = contract(T, [1,4,3], EnV.A{i+1}, [4,1,2]);
        elseif i == length(EnV.A) - 1
            T = contract(T, 1, EnV.A{i-1}, 3);
            T = contract(T, 2, EnV.A{i+1}, 3);
            T = contract(T, [1,3,4], EnV.A{i}, [4,1,2], [1,3,2]);
        else
            T = contract(T, 1, EnV.A{i-1}, 3);
            T = contract(T, [2,5], EnV.A{i}, [4,1]);
            T = contract(T, [2,4], EnV.A{i+1}, [4,1]);
            T = contract(T, [1,4], EnV.A{i+2}, [3,2], [1,4,2,3]);
        end
    case '-<-'
        if i == length(EnV.A)
            T = contract(T, 1, EnV.A{i-2}, 3);
            T = contract(T, 2, EnV.A{i}, 3);
            T = contract(T, [1,3,4], EnV.A{i-1}, [4,1,2], [1,3,2]);
        elseif i == 2
            T = contract(T, 1, EnV.A{i+1}, 3);
            T = contract(T, 1, EnV.A{i-1}, 3);
            T = contract(T, [1,4,3], EnV.A{i}, [4,1,2]);
        else
            T = contract(T, 1, EnV.A{i-2}, 3);
            T = contract(T, [2,5], EnV.A{i-1}, [4,1]);
            T = contract(T, [2,4], EnV.A{i}, [4,1]);
            T = contract(T, [1,4], EnV.A{i+1}, [3,2], [1,4,2,3]);
        end
end
end

