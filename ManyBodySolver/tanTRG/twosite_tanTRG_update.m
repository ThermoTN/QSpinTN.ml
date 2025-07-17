function [ rho, EnV, TE, EE ] = twosite_tanTRG_update(EnV, rho, pos, dir, tau, Para)
D_max = Para.MCrit;
switch dir
    case '->-'
        if pos == 1
            T = contract(rho.A{pos}, 1, rho.A{pos+1}, 1, [3,1,4,2,5]);
        elseif pos == length(EnV.A) - 1
            T = contract(rho.A{pos}, 2, rho.A{pos+1}, 1, [1,2,4,3,5]);
        else
            T = contract(rho.A{pos}, 2, rho.A{pos+1}, 1, [1,4,2,5,3,6]);
        end
    case '-<-'
        if pos == length(EnV.A)
            T = contract(rho.A{pos-1}, 2, rho.A{pos}, 1, [1,2,4,3,5]);
        elseif pos == 2
            T = contract(rho.A{pos-1}, 1, rho.A{pos}, 1, [3,1,4,2,5]);
        else
            T = contract(rho.A{pos-1}, 2, rho.A{pos}, 1, [1,4,2,5,3,6]);
        end
end

fun_prod = @(T) twosite_prod(T, EnV, pos, dir);
fun_inner = @(T1, T2) tensor_inner_prod(T1, T2);
fun_norm = @(T) tensor_norm(T);
[T, Info] = LanczosExp(fun_prod, T, -tau, 'innerprod', fun_inner, 'norm', fun_norm);


len = length(rho.A);
if pos == 1 || (pos == 2 && strcmp(dir, '-<-'))
    T = permute(T, [2,4,1,3,5]);
    Tsize = size(T);
    T = reshape(T, [prod(Tsize(1:1:2)), prod(Tsize(3:1:5))]);
elseif (pos == (len-1) && strcmp(dir, '->-')) || pos == len
    T = permute(T, [1,2,4,3,5]);
    Tsize = size(T);
    T = reshape(T, [prod(Tsize(1:1:3)), prod(Tsize(4:1:5))]);
else
    T = permute(T, [1,3,5,2,4,6]);
    Tsize = size(T);
    T = reshape(T, [prod(Tsize(1:1:3)), prod(Tsize(4:1:6))]);
end

[U, S, V, TE, EE] = svdT(T, 'Nkeep', D_max, 'epsilon', Para.TSRGStol);
Usize = size(U);
Vsize = size(V);
if pos == 1 || (pos == 2 && strcmp(dir, '-<-'))
    U = reshape(U, [Tsize(1), Tsize(2), Usize(2)]);
    U = permute(U, [3,1,2]);
    V = reshape(V, [Vsize(1), Tsize(3), Tsize(4), Tsize(5)]);
elseif (pos == (len-1) && strcmp(dir, '->-')) || pos == len
    U = reshape(U, [Tsize(1), Tsize(2), Tsize(3), Usize(2)]);
    U = permute(U, [1,4,2,3]);
    V = reshape(V, [Vsize(1), Tsize(4), Tsize(5)]);
else
    U = reshape(U, [Tsize(1), Tsize(2), Tsize(3), Usize(2)]);
    U = permute(U, [1,4,2,3]);
    V = reshape(V, [Vsize(1), Tsize(4), Tsize(5), Tsize(6)]);
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
            rho.A{pos-1} = contract(U, 1, S, 1, [3,1,2]);
        else
            rho.A{pos-1} = contract(U, 2, S, 1, [1,4,2,3]);
        end
        rho.A{pos} = V;
end
log(Ns);
rho.lgnorm = rho.lgnorm + log(Ns);

% updata EnV
switch dir
    case '->-'
        if pos == 1
            EnV.A{pos} = contract(rho.A{pos}, 2, EnV.A{pos}, 3);
            EnV.A{pos} = contract(EnV.A{pos}, [2,4], conj(rho.A{pos}), [3,2], [3,2,1]);
        else
            EnV.A{pos} = contract(EnV.A{pos-1}, 2, EnV.A{pos}, 1);
            EnV.A{pos} = contract(EnV.A{pos}, [1,4], conj(rho.A{pos}), [1,3]);
            EnV.A{pos} = contract(EnV.A{pos}, [1,3,5], rho.A{pos}, [1,3,4], [2,1,3]);
        end
    case '-<-'
        if pos == length(EnV.A)
            EnV.A{pos} = contract(rho.A{pos}, 2, EnV.A{pos}, 3);
            EnV.A{pos} = contract(EnV.A{pos}, [2,4], conj(rho.A{pos}), [3,2], [3,2,1]);
        else
            EnV.A{pos} = contract(EnV.A{pos+1}, 2, EnV.A{pos}, 2);
            EnV.A{pos} = contract(EnV.A{pos}, [1,4], conj(rho.A{pos}), [2,3]);
            EnV.A{pos} = contract(EnV.A{pos}, [1,3,5], rho.A{pos}, [2,3,4], [2,1,3]);
        end
end
end

function [ T ] = twosite_prod(T, EnV, i, dir)
switch dir
    case '->-'
        if i == 1
            T = contract(T, 1, EnV.A{i+2}, 3);
            T = contract(T, 1, EnV.A{i}, 3);
            T = contract(T, [1,6,5], EnV.A{i+1}, [4,1,2], [3,4,5,1,2]);
        elseif i == length(EnV.A) - 1
            T = contract(T, 1, EnV.A{i-1}, 3);
            T = contract(T, 2, EnV.A{i+1}, 3);
            T = contract(T, [1,5,6], EnV.A{i}, [4,1,2], [3,5,4,1,2]);
        else
            T = contract(T, 1, EnV.A{i-1}, 3);
            T = contract(T, [2,7], EnV.A{i}, [4,1]);
            T = contract(T, [2,6], EnV.A{i+1}, [4,1]);
            T = contract(T, [1,6], EnV.A{i+2}, [3,2], [3,6,4,5,1,2]);
        end
    case '-<-'
        if i == length(EnV.A)
            T = contract(T, 1, EnV.A{i-2}, 3);
            T = contract(T, 2, EnV.A{i}, 3);
            T = contract(T, [1,5,6], EnV.A{i-1}, [4,1,2], [3,5,4,1,2]);
        elseif i == 2
            T = contract(T, 1, EnV.A{i+1}, 3);
            T = contract(T, 1, EnV.A{i-1}, 3);
            T = contract(T, [1,6,5], EnV.A{i}, [4,1,2], [3,4,5,1,2]);
        else
            T = contract(T, 1, EnV.A{i-2}, 3);
            T = contract(T, [2,7], EnV.A{i-1}, [4,1]);
            T = contract(T, [2,6], EnV.A{i}, [4,1]);
            T = contract(T, [1,6], EnV.A{i+1}, [3,2], [3,6,4,5,1,2]);
        end
end
end

