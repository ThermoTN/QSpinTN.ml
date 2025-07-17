function [ EnV ] = Init_tdvp_EnV(H, T)
EnV = H;
Vr = H.A{end};

for i = (length(H.A)-1):-1:2
    if i == (length(H.A)-1)
        Vr = contract(T.A{i+1}, 2, Vr, 3);
        Vr = contract(conj(T.A{i+1}), 2, Vr, 3);
        EnV.A{i+1} = permute(Vr, [1,3,2]);
        Vr = contract(Vr, 3, H.A{i},2);
    else
        Vr = contract(T.A{i+1}, [2,3], Vr, [2,5]);
        Vr = contract(conj(T.A{i+1}), [2,3], Vr, [2,4]);
        EnV.A{i+1} = permute(Vr, [1,3,2]);
        Vr = contract(Vr, 3, H.A{i}, 2);
    end
end
end