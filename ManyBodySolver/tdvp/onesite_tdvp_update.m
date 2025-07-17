function [ rho, EnV, R ] = onesite_tdvp_update(EnV, rho, pos, dir, tau, mode)
T = rho.A{pos};
fun_prod = @(T) onesite_prod(T, EnV, pos);
fun_inner = @(T1, T2) tensor_inner_prod(T1, T2);
fun_norm = @(T) tensor_norm(T);
[T, Info] = LanczosExp(fun_prod, T, tau, 'innerprod', fun_inner, 'norm', fun_norm);
% [T, Info] = LanczosExp(fun_prod, T, -tau, 'innerprod', fun_inner, 'norm', fun_norm);
switch dir
    case '->-'
        if pos == 1
            sizeT = size(T);
            T = permute(T, [2,1]);
            T = reshape(T, [sizeT(2), sizeT(1)]);
            [Q, S, V] = svd(T, 'econ');
            rho.A{pos} = permute(Q, [2,1]);
            Ns = norm(diag(S));
            R = S * V';
        elseif pos == length(EnV.A)
            rho.A{pos} = T;
            sizeT = size(T);
            T = permute(T, [2,1]);
            T = reshape(T, [sizeT(2), sizeT(1)]);
            [~, R, ~] = svd(T, 'econ');
            Ns = norm(diag(R));
            rho.A{pos} = rho.A{pos} / Ns;
        else
            sizeT = size(T);
            T = permute(T, [1,3,2]);
            T = reshape(T, [sizeT(1) * sizeT(3), sizeT(2)]);
            [Q, S, V] = svd(T, 'econ');
            Q = reshape(Q, [sizeT(1), sizeT(3), size(S, 1)]);
            rho.A{pos} = permute(Q, [1,3,2]);
            Ns = norm(diag(S));
            R = S * V';
        end
    case '-<-'
        if pos == length(EnV.A)
            %             sizeT = size(T);
            %             T = permute(T, [2,3,1]);
            %             T = reshape(T, [sizeT(2) * sizeT(3), sizeT(1)]);
            %             [Q, S, V] = svd(T, 'econ');
            %             Q = reshape(Q, [sizeT(2), sizeT(3), sizeT(1)]);
            %             rho.A{pos} = permute(Q, [3,1,2]);
            %             Ns = norm(diag(S));
            %             R = S * V';
            %             R = permute(R, [2,1]);
            sizeT = size(T);
            [U, S, V] = svd(T, 'econ');
            rho.A{pos} = V';
            Ns = norm(diag(S));
            R = U * S;
        elseif pos == 1
            rho.A{pos} = T;
            sizeT = size(T);
            T = permute(T, [2,1]);
            T = reshape(T, [sizeT(2), sizeT(1)]);
            [~, R, ~] = svd(T, 'econ');
            Ns = norm(diag(R));
            rho.A{pos} = rho.A{pos} / Ns;
            % Ns = log(norm(R));
        else
            %             sizeT = size(T);
            %             T = permute(T, [2,3,4,1]);
            %             T = reshape(T, [sizeT(2) * sizeT(3) * sizeT(4), sizeT(1)]);
            %             [Q, S, V] = svd(T, 'econ');
            %             Q = reshape(Q, [sizeT(2), sizeT(3), sizeT(4), sizeT(1)]);
            %             rho.A{pos} = permute(Q, [4,1,2,3]);
            %             Ns = norm(diag(S));
            %             R = S * V';
            %             R = permute(R, [2,1]);
            sizeT = size(T);
            T = reshape(T, [sizeT(1), sizeT(2) * sizeT(3)]);
            [U, S, V] = svd(T, 'econ');
            rho.A{pos} = reshape(V', [size(S, 1), sizeT(2), sizeT(3)]);
            Ns = norm(diag(S));
            R = U * S;
        end
end
R = R/Ns;
log(Ns);
rho.lgnorm = rho.lgnorm + log(Ns);


switch mode
    case '2t1'
        % update rho
        switch dir
            case '->-'
                if pos == 1
                    error('1!')
                elseif pos == length(rho.A)
                    % do nothing
                else
                    rho.A{pos+1} = contract(R, 2, rho.A{pos+1}, 1);
                end
            case '-<-'
                if pos == length(rho.A)
                    error('1!')
                elseif pos == 2
                    rho.A{pos-1} = contract(R, 1, rho.A{pos-1}, 1);
                elseif pos == 1
                    % do nothing
                else
                    rho.A{pos-1} = contract(R, 1, rho.A{pos-1}, 2, [2,1,3]);
                end
        end
        R = 0;
        
    case '1t0'
        % update EnV
        switch dir
            case '->-'
                if pos == 1
                    EnV.A{pos} = contract(rho.A{pos}, 2, EnV.A{pos}, 3);
                    EnV.A{pos} = contract(EnV.A{pos}, 3, conj(rho.A{pos}), 2, [3,2,1]);
                elseif pos < length(EnV.A)
                    EnV.A{pos} = contract(EnV.A{pos-1}, 2, EnV.A{pos}, 1);
                    EnV.A{pos} = contract(EnV.A{pos}, [1,4], conj(rho.A{pos}), [1,3]);
                    EnV.A{pos} = contract(EnV.A{pos}, [1,3], rho.A{pos}, [1,3], [2,1,3]);
                end
            case '-<-'
                if pos == length(EnV.A)
                    EnV.A{pos} = contract(rho.A{pos}, 2, EnV.A{pos}, 3);
                    EnV.A{pos} = contract(EnV.A{pos}, 3, conj(rho.A{pos}), 2, [3,2,1]);
                elseif pos > 1
                    EnV.A{pos} = contract(EnV.A{pos+1}, 2, EnV.A{pos}, 2);
                    EnV.A{pos} = contract(EnV.A{pos}, [1,4], conj(rho.A{pos}), [2,3]);
                    EnV.A{pos} = contract(EnV.A{pos}, [1,3], rho.A{pos}, [2,3], [2,1,3]);
                end
        end
end
end

function [ T ] = onesite_prod(T, EnV, i)

if i == 1
    T = contract(T, 1, EnV.A{i+1}, 3);
    T = contract(T, [1,3], EnV.A{i}, [3,1]);
elseif i == length(EnV.A)
    T = contract(T, 1, EnV.A{i-1}, 3);
    T = contract(T, [1,3], EnV.A{i}, [3,1]);
else
    T = contract(T, 1, EnV.A{i-1}, 3);
    T = contract(T, [2,4], EnV.A{i}, [4,1]);
    T = contract(T, [1,3], EnV.A{i+1}, [3,2], [1,3,2]);
end

end

