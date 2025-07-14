function [ n ] = tensor_norm(T)
n = norm(reshape(T, [numel(T), 1]));
end