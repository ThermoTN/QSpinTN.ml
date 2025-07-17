function [A] = compressMPS(A)
L = length(A);
for it = L:-1:2
    if it == L
        [U, S, V] = svd(A{it}, 'econ');
        A{it} = V';
        R = U * S;
        A{it-1} = contract(A{it-1}, 2, R, 1, [1,3,2]);
    else
        sizeA = size(A{it});
        T = reshape(A{it}, [sizeA(1), sizeA(2) * sizeA(3)]);
        [U, S, V] = svd(T, 'econ');
        A{it} = reshape(V', [size(S,1), sizeA(2), sizeA(3)]);
        R = U * S;
        if it == 2
            A{it-1} = contract(A{it-1}, 1, R, 1, [2,1]);
        else
            A{it-1} = contract(A{it-1}, 2, R, 1, [1,3,2]);
        end
    end
    
end
end

