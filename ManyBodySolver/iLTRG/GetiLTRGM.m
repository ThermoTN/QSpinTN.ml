function M = GetiLTRGM( Ta, La, Tb, Lb, Para )

% Get left environment
Tl = eye(size(La));
fmaxval = 0;

for It = 1:1:1000
    Tl = TransM(Tl, Ta, La, 'L2R');
    Tl = TransM(Tl, Tb, Lb, 'L2R');
    
    maxval = norm(Tl);
    Tl = Tl/maxval;
    
    if abs(maxval-fmaxval) < 1e-9
        break;
    end
    fmaxval = maxval;
end

Tr = eye(size(La));
fmaxval = 0;

for It = 1:1:1000
    Tr = TransM(Tr, Tb, Lb, 'R2L');
    Tr = TransM(Tr, Ta, La, 'R2L');
    
    maxval = norm(Tr);
    Tr = Tr/maxval;
    
    if abs(maxval-fmaxval) < 1e-9
        break;
    end
    fmaxval = maxval;
end

h = Para.Field.h/norm(Para.Field.h);
[Sx, Sy, Sz, Id] = SpinOp(Para.d);

MOp = h(1) * Sx + h(2) * Sy + h(3) * Sz;

Tap = contract(Ta, 3, MOp, 2, [1,2,4,3]);

Tmp = contract(Tl, 1, La, 1);
Tmp = contract(Tmp, 1, conj(La), 1);
Tmp = contract(Tmp, 1, Tap, 1);
Tmp = contract(Tmp, [1,3,4], conj(Ta), [1,3,4]);
Tmp = TransM(Tmp, Tb, Lb, 'L2R');
Ma = contract(Tmp, [1,2], Tr, [1,2]);

Tbp = contract(Tb, 3, MOp, 2, [1,2,4,3]);
Tmp = TransM(Tl, Ta, La, 'L2R');

Tmp = contract(Tmp, 1, Lb, 1);
Tmp = contract(Tmp, 1, conj(Lb), 1);
Tmp = contract(Tmp, 1, Tbp, 1);
Tmp = contract(Tmp, [1,3,4], conj(Tb), [1,3,4]);
Mb = contract(Tmp, [1,2], Tr, [1,2]);
M = (Ma + Mb)/(2 * contract(Tl, [1,2], Tr, [1,2]));
end

