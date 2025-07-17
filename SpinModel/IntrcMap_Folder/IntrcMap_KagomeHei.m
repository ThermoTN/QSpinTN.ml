function [Intr] = IntrcMap_KagomeHei(Para)

Lx = Para.Geo.Lx;
Ly = Para.Geo.Ly;
BCY = Para.Geo.BCY;
int_num = 6 * Lx * Ly - 2 * (Lx + Ly) + 1;
if strcmp(BCY, 'PBC')
    int_num = int_num + 2 * Lx - 1;
end
int_cell = cell(int_num, 3);
int_count = 1;
for Lx_i = 1:1:Lx
    base_num = 3 * Ly * (Lx_i - 1);
    
    for i = 1:1:(2 * Ly - 1)
        int_cell{int_count, 1} = base_num + i;
        int_cell{int_count, 2} = base_num + i + 1;
        int_cell{int_count, 3} = 'Hei';
        int_count = int_count + 1;
    end
    
    for i = 1:1:Ly
        int_cell{int_count, 1} = base_num + 2 * i - 1;
        int_cell{int_count, 2} = base_num + 2 * Ly + i;
        int_cell{int_count, 3} = 'Hei';
        int_count = int_count + 1;
        
        int_cell{int_count, 1} = base_num + 2 * i;
        int_cell{int_count, 2} = base_num + 2 * Ly + i;
        int_cell{int_count, 3} = 'Hei';
        int_count = int_count + 1;
    end
    
    if Lx_i < Lx
        for i = 1:1:Ly
            int_cell{int_count, 1} = base_num + 2 * Ly + i;
            int_cell{int_count, 2} = base_num + 3 * Ly + 2 * i;
            int_cell{int_count, 3} = 'Hei';
            int_count = int_count + 1;
            
            if i < Ly
                int_cell{int_count, 1} = base_num + 2 * Ly + i;
                int_cell{int_count, 2} = base_num + 3 * Ly + 2 * i + 1;
                int_cell{int_count, 3} = 'Hei';
                int_count = int_count + 1;
            end
        end
    end
    
    if strcmp(BCY, 'PBC')
        int_cell{int_count, 1} = base_num + 1;
        int_cell{int_count, 2} = base_num + 2 * Ly;
        int_cell{int_count, 3} = 'Hei';
        int_count = int_count + 1;
        if Lx_i < Lx
            int_cell{int_count, 1} = base_num + 3 * Ly;
            int_cell{int_count, 2} = base_num + 3 * Ly + 1;
            int_cell{int_count, 3} = 'Hei';
            int_count = int_count + 1;
        end
    end
end
int_count - 1;
int_cell = GetIntCell(int_cell, Para);
Intr = struct('JmpOut', int_cell(:,1), 'JmpIn', int_cell(:,2), ...
              'JmpOut_type', int_cell(:,3), 'JmpIn_type', int_cell(:,4), 'CS', int_cell(:,5));

for i = 1:1:length(Intr)
    if Intr(i).JmpOut > Intr(i).JmpIn
        T = Intr(i).JmpIn;
        Intr(i).JmpIn = Intr(i).JmpOut;
        Intr(i).JmpOut = T;
    end
end
% keyboard;
end

function int_cell = GetIntCell(old_int_cell, Para)
int_cell = cell(3 * length(old_int_cell(:,1)), 5);
spin_cell = {'Sy', 'Sz', 'Sx'};
for i = 1:length(old_int_cell(:,1))
    for j = 3*(i-1)+1:3*i
        int_cell(j, 1) = old_int_cell(i,1);
        int_cell(j, 2) = old_int_cell(i,2);
        int_cell{j, 3} = spin_cell{1 + mod(j,3)};
        int_cell{j, 4} = spin_cell{1 + mod(j,3)};
        int_cell{j, 5} = Para.Model.J;
    end
end
end