function [ Bond_cell, nn_num, nnn_num ] = Bond_Info_Kitaev( Para )

Lx = Para.Geo.Lx;
Ly = Para.Geo.Ly;
% keyboard;
if mod(Lx, 2) ~= 1 || mod(Ly, 2) ~= 0
    error('Wrong geometry!')
end
L = Lx * Ly;
Bond_num = (Lx - 1) * (Ly - 1) / 2 + Ly / 2 + Ly * (Lx - 1);
BCY = Para.Geo.BCY;


if strcmp(BCY, 'PBC')
    Bond_num = Bond_num + (Lx - 1) / 2;
end
Bond_cell = cell(Bond_num, 3);
count = 1;
nn_num = 0;
nnn_num = 0;
% OBC part

    for x = 1:1:Lx
    bnum = Ly * (x - 1);
    % z-type
    if mod(x, 2) == 1
        for y = 1:2:(Ly-1)
            Bond_cell{count, 1} = bnum + y;
            Bond_cell{count, 2} = bnum + y + 1;
            Bond_cell{count, 3} = 'z-type';
            
            count = count + 1;
        end
    else
        for y = 2:2:(Ly-2)
            Bond_cell{count, 1} = bnum + y;
            Bond_cell{count, 2} = bnum + y + 1;
            Bond_cell{count, 3} = 'z-type';
            
            count = count + 1;
        end
    end
    
    % x-type
    if mod(x, 2) == 1 && x ~= Lx
        for y = 1:1:Ly
            Bond_cell{count, 1} = bnum + y;
            Bond_cell{count, 2} = bnum + y + Ly;
            Bond_cell{count, 3} = 'x-type';
            
            count = count + 1;
        end
    end
    
    % y-type
    if mod(x,2) == 0
        for y = 1:1:Ly
            Bond_cell{count, 1} = bnum + y;
            Bond_cell{count, 2} = bnum + y + Ly;
            Bond_cell{count, 3} = 'y-type';
            
            count = count + 1;
        end
    end
end

if strcmp(BCY, 'PBC')
    for x = 2:2:Lx
        Bond_cell{count, 1} = Ly * (x - 1) + 1;
        Bond_cell{count, 2} = Ly * x;
        Bond_cell{count, 3} = 'z-type';
        
        count = count + 1;
    end
    
end

end