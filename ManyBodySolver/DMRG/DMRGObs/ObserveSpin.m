function [ObSpin] = ObserveSpin(S, GS)
% Observation of <Sa>
d = 2;
L = length(GS);
ObSpin = zeros(L, 1);
left = cell(L-1, 1);
%  O---2
%  |
%  O---1

for  i = 1:L
    if i == 1
        ObSpini = contract(S, 1, GS{i}, 2);
        ObSpini = contract(conj(GS{i}), 2, ObSpini, 1);
        bond1 = length(ObSpini(:,1));
        bond2 = length(ObSpini(1,:));
        ObSpini = contract(ObSpini, [1,2], eye(bond1, bond2), [1,2]);
        ObSpin(i) = ObSpini;
        left{i} = contract(GS{1}, 2, conj(GS{1}), 2);
        
    else
        if i ~= L
            left{i} = contract(left{i-1}, 2, conj(GS{i}), 1);
            left{i} = contract(GS{i}, 1, left{i}, 1);
            ObSpini = contract(left{i}, [2,4], S, [1,2]);
            bond = length(GS{i}(1,:,1));
            ObSpin(i) = contract(ObSpini, [1,2], eye(bond, bond), [1,2]);
            left{i} = contract(left{i}, [2,4], eye(d, d), [1,2]);
        else
            left{i} = contract(left{i-1}, 2, conj(GS{i}), 1);
            left{i} = contract(GS{i}, 1, left{i}, 1);
            ObSpin(i) = contract(left{i}, [1,2], S, [1,2]);
        end
    end
end

end