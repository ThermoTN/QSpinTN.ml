function T = GetEn(EnV, H, rho)
T = contract(EnV, 1, conj(rho.A{1}), 1);
T = contract(T, [1,3], H.A{1}, [1,2]);
T = contract(T, [1,3,2], rho.A{1}, [1,2,3]);
end