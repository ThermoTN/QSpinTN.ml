function [ H ] = InitSM(PMPO, Para)

H = PMPO(1).A;
for it = 2:Para.L
    H = kronS(H, PMPO(it).A);
end
H = H{1};
end



function C = kronS(A, B)

szA = size(A);
szB = size(B);

C = cell([szA(1), szB(2)]);

sz0 = size(kron(A{1, 1}, B{1, 1}));
for iti = 1:szA(1)
    for itj = 1:szB(2)
        C{iti, itj} = sparse(sz0(1), sz0(2));
        for its = 1:szA(2)
            C{iti, itj} = C{iti, itj} + kron(A{iti, its}, B{its, itj});
        end
    end
end

end