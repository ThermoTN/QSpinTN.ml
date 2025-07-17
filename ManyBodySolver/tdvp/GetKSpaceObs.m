function [ Obs ] = GetKSpaceObs(Qrho, rho, QuOp)

QuOp.A = GetDagger(QuOp.A);
L = length(Qrho.A);
Vr = QuOp.A{end};

for i = (L-1):-1:1
    if i == (L-1)
        Vr = contract(Qrho.A{i+1}, 2, Vr, 3);
        Vr = contract(conj(rho.A{i+1}), 2, Vr, 3);
        Vr = contract(Vr, 3, QuOp.A{i},2);
    elseif i == 1
        Vr = contract(Qrho.A{i+1}, [2,3], Vr, [2,5]);
        Vr = contract(conj(rho.A{i+1}), [2,3], Vr, [2,4]);
        Vr = contract(Vr, 3, QuOp.A{i}, 1);
    else
        Vr = contract(Qrho.A{i+1}, [2,3], Vr, [2,5]);
        Vr = contract(conj(rho.A{i+1}), [2,3], Vr, [2,4]);
        Vr = contract(Vr, 3, QuOp.A{i}, 2);
    end
end
Vr = contract(Qrho.A{1}, [1,2], Vr, [2,4]);
Vr = contract(conj(rho.A{1}), [1,2], Vr, [1,2]);
Obs = Vr;
end

function A = GetDagger(A)
for it = 1:length(A)
    if it == 1 || it == length(A)
        A{it} = permute(conj(A{it}), [1,3,2]);
    else
        A{it} = permute(conj(A{it}), [1,2,4,3]);
    end
end
end