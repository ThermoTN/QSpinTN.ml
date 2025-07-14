function En = GetEn(EnV, T)
Vl = EnV.A{1};
Vr = EnV.A{2};
TT = contract(T.A{1}, 1, T.A{2}, 1);
TTp = contract(TT, 1, Vl, 3);
TTp = contract(TTp, [1,2,3], Vr, [2,5,3]);
En = contract(TT,[1,2,3], conj(TTp), [1,2,3]);
end