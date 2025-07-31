function [ T ] = TransM( T, Ta, La, ChLR )

switch ChLR
    case 'L2R'
        T = contract(T, 1, La, 1);
        T = contract(T, 1, conj(La), 1);
        T = contract(T, 1, Ta, 1);
        T = contract(T, [1,3,4], conj(Ta), [1,3,4]);
    case 'R2L'
        T = contract(T, 1, Ta, 2);
        T = contract(T, [1,3,4], conj(Ta), [2,3,4]);
        T = contract(T, 1, La, 2);
        T = contract(T, 1, conj(La), 2);
end
end

