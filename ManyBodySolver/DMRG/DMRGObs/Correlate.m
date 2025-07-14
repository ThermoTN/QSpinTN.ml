function [CorSpin] = Correlate(Si, Sj, GS)
% 2022/4/14 by ezlv
% Corralate <SaSb> Sa is on i and Sb is on j
% a, b = x, y, z
L = length(GS);
d = 2;
CorSpin = zeros(L, L);
% MPS{i} = contract i tensors of GS
%  O---2
%  |
%  O---1

for i = 1:L % Si
    for j = i+1:L % Sj
	left = cell(L-1, 1);
        for k = 1:j
            if i == 1 && k == 1
		% keyboard;
                left{k} = contract(GS{k}, 2, Si, 1);
                left{k} = contract(left{k}, 2, conj(GS{k}), 2);
            elseif j == L && k == L
		% keyboard;
                left{k} = contract(GS{k}, 1, left{k-1}, 1);
                left{k} = contract(left{k}, 2, conj(GS{k}), 1);
                CorSpin(i,j) = contract(left{k}, [1,2], eye(d,d), [1,2]);
            else
                if k == 1
                    left{k} = contract(GS{k}, 2, conj(GS{k}), 2);
                elseif k == i
                    left{k} = contract(GS{k}, 1, left{k-1}, 1);
                    left{k} = contract(left{k}, 3, conj(GS{k}), 1);
                    left{k} = contract(left{k}, [2,4], Si, [1,2]);
                elseif k == j
                    % keyboard;
                    left{k} = contract(GS{k}, 1, left{k-1}, 1);
                    left{k} = contract(left{k}, 3, conj(GS{k}), 1);
                    left{k} = contract(left{k}, [2,4], Sj, [1,2]);
                    bond = length(GS{k}(1,:,1));
                    CorSpin(i,j) = contract(left{k}, [1,2], eye(bond, bond), [1,2]);
                else
                    left{k} = contract(GS{k}, 1, left{k-1}, 1);
                    left{k} = contract(left{k}, 3, conj(GS{k}), 1);
                    left{k} = contract(left{k}, [2,4], eye(d,d), [1,2]);
                end
            end 
        end
    end
end
end
