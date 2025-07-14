function [Intr] = IntrcMap_ToricCode(Para)
% L = 4;
L = Para.Geo.LT;
d = Para.d;
Lx = 2*L; Ly = L;

[Sx, ~, Sz, ~] = SpinOp(d);



cnt = 0;
A_term = [];
for x = 3:2:Lx+1
    for y = 1:Ly
        pos1 = get_idx_snake_PBC(x, y, Lx, Ly);
        pos2 = get_idx_snake_PBC(x+1, y+1, Lx, Ly);
        pos3 = get_idx_snake_PBC(x, y+1, Lx, Ly);
        pos4 = get_idx_snake_PBC(x-1, y+1, Lx, Ly);
        %         A_term = [A_term; sort([pos1, pos2, pos3, pos4])];
        cnt = cnt + 1;
        Intr(cnt).site = sort([pos1, pos2, pos3, pos4]);
        Intr(cnt).J = Para.Model.JA;
        Intr(cnt).operator = {Sx, Sx, Sx, Sx};
    end
end




B_term = [];
for x = 2:2:Lx
    for y = 1:Ly
        pos1 = get_idx_snake_PBC(x, y, Lx, Ly);
        pos2 = get_idx_snake_PBC(x+1, y, Lx, Ly);
        pos3 = get_idx_snake_PBC(x, y+1, Lx, Ly);
        pos4 = get_idx_snake_PBC(x-1, y, Lx, Ly);
        %         B_term = [B_term; sort([pos1, pos2, pos3, pos4])];
        cnt = cnt + 1;
        Intr(cnt).site = sort([pos1, pos2, pos3, pos4]);
        Intr(cnt).J = Para.Model.JB;
        Intr(cnt).operator = {Sz, Sz, Sz, Sz};
    end
end

for it = 1:length(Intr)
    OO_list(it) = Intr(it).site(1);
end

[~, ord] = sort(OO_list);

Intr = Intr(ord);
end

function idx = get_idx_snake_PBC(x, y, Lx, Ly)
if x > Lx
    x = mod(x, Lx);
end
if y > Ly
    y = mod(y, Ly);
end
if mod(x, 2) == 1 % Odd columns (x = 1, 3, ...), y goes from 1 to Ly
    idx = (x - 1) * Ly + y;
else % Even columns (x = 2, 4, ...), y goes from Ly down to 1
    idx = (x - 1) * Ly + (Ly - y + 1);
end
end