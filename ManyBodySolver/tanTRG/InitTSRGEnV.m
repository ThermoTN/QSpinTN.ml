function [ EnV ] = InitTSRGEnV(rho, H)
EnV = H;
EnV.A{end} = contract(H.A{end}, 2, conj(rho.A{end}), 2, [3,1,4,2]);
EnV.A{end} = contract(EnV.A{end}, [3,4], rho.A{end}, [3,2]);

for i = (length(H.A)-1):-1:2
    EnV.A{i} = contract(EnV.A{i+1}, 1, conj(rho.A{i}), 2);
    EnV.A{i} = contract(EnV.A{i}, [1,4], H.A{i}, [2,3]);
    EnV.A{i} = contract(EnV.A{i}, [1,5,3], rho.A{i}, [2,3,4]);
end
end