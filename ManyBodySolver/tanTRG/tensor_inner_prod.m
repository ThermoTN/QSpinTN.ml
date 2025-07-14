function [ n ] = tensor_inner_prod(T1, T2)
n = reshape(T1, [numel(T1), 1])' * reshape(T2, [numel(T2), 1]);
end