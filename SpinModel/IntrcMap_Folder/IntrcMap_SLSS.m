function [ Intr ] = IntrcMap_SLSS( Para )
% Triangular lattice 
% Transverse field Ising model
% Parameter: 
%           J1     Nearest neighbor term
%           J2     Next-nearest neighbor term
%           Delta  Transverse field term
% 
% Hamiltonian:
%   H = J1 \sum_<i,j> Sz_i Sz_j 
%       + J2 \sum_<<i,j>> Sz_i Sz_j 
%       - Delta\sum_i Sx_i 
%       - hz\sum_i Sz_i
Lx = Para.Geo.Lx;
Ly = Para.Geo.Ly;
L = Lx * Ly;
int_num = Lx*(Ly-1)+Ly*(Lx-1)+floor(Lx/2)*floor(Ly/2) + floor((Lx-1)/2) * floor((Ly-1)/2);
BCX = Para.Geo.BCX;
BCY = Para.Geo.BCY;
if strcmp(BCX, 'PBC')
    int_num = int_num + Ly + floor(Ly/2)*(floor((Lx+1)/2)-floor(Lx/2)) + floor((Ly-1)/2)*(floor(Lx/2)-floor((Lx-1)/2));
end

if strcmp(BCY, 'PBC')
    int_num = int_num + Lx + floor(Lx/2)*(floor((Ly+1)/2)-floor(Ly/2)) + floor((Lx-1)/2)*(floor(Ly/2)-floor((Ly-1)/2));
end
int_num = 3 * int_num;
int_cell = cell(int_num, 5);
count = 1;

% OBC part

for x = 1:1:Lx
    bnum = Ly * (x - 1);
    
    % J along y dir
    for y = 1:1:(Ly-1)
        
        int_cell{count, 1} = bnum + y;
        int_cell{count, 2} = bnum + y + 1;
        int_cell{count, 3} = 'Sx';
        int_cell{count, 4} = 'Sx';
        int_cell{count, 5} = Para.Model.J;
        
        int_cell{count+1, 1} = bnum + y;
        int_cell{count+1, 2} = bnum + y + 1;
        int_cell{count+1, 3} = 'Sy';
        int_cell{count+1, 4} = 'Sy';
        int_cell{count+1, 5} = Para.Model.J;
        
        int_cell{count+2, 1} = bnum + y;
        int_cell{count+2, 2} = bnum + y + 1;
        int_cell{count+2, 3} = 'Sz';
        int_cell{count+2, 4} = 'Sz';
        int_cell{count+2, 5} = Para.Model.J;
        
        count = count + 3;
    end
    
    if x ~= Lx
        % J along x dir
        for y = 1:1:Ly
        
            int_cell{count, 1} = bnum + y;
            int_cell{count, 2} = bnum + Ly + y;
            int_cell{count, 3} = 'Sx';
            int_cell{count, 4} = 'Sx';
            int_cell{count, 5} = Para.Model.J;
            
            int_cell{count+1, 1} = bnum + y;
            int_cell{count+1, 2} = bnum + Ly + y;
            int_cell{count+1, 3} = 'Sy';
            int_cell{count+1, 4} = 'Sy';
            int_cell{count+1, 5} = Para.Model.J;
            
            int_cell{count+2, 1} = bnum + y;
            int_cell{count+2, 2} = bnum + Ly + y;
            int_cell{count+2, 3} = 'Sz';
            int_cell{count+2, 4} = 'Sz';
            int_cell{count+2, 5} = Para.Model.J;
            
            count = count + 3;
        end
        
    end
    
    for y = 1:1:(Ly-1)
        % Jp: '\'
        if mod(x,2) == 1 && mod(y,2) == 1 && x ~= Lx
            int_cell{count, 1} = bnum + y;
            int_cell{count, 2} = bnum + Ly + y + 1;
            int_cell{count, 3} = 'Sx';
            int_cell{count, 4} = 'Sx';
            int_cell{count, 5} = Para.Model.JD;
            
            int_cell{count+1, 1} = bnum + y;
            int_cell{count+1, 2} = bnum + Ly + y + 1;
            int_cell{count+1, 3} = 'Sy';
            int_cell{count+1, 4} = 'Sy';
            int_cell{count+1, 5} = Para.Model.JD;
            
            int_cell{count+2, 1} = bnum + y;
            int_cell{count+2, 2} = bnum + Ly + y + 1;
            int_cell{count+2, 3} = 'Sz';
            int_cell{count+2, 4} = 'Sz';
            int_cell{count+2, 5} = Para.Model.JD;
            
            count = count + 3;
        end
    end
    
    for y = 3:1:(Ly)
        % Jp: '/'
        if mod(x,2) == 0 && mod(y,2) == 1 && x ~= Lx
            int_cell{count, 1} = bnum + y;
            int_cell{count, 2} = bnum + Ly + y - 1;
            int_cell{count, 3} = 'Sx';
            int_cell{count, 4} = 'Sx';
            int_cell{count, 5} = Para.Model.JD;
            
            int_cell{count+1, 1} = bnum + y;
            int_cell{count+1, 2} = bnum + Ly + y - 1;
            int_cell{count+1, 3} = 'Sy';
            int_cell{count+1, 4} = 'Sy';
            int_cell{count+1, 5} = Para.Model.JD;
            
            int_cell{count+2, 1} = bnum + y;
            int_cell{count+2, 2} = bnum + Ly + y - 1;
            int_cell{count+2, 3} = 'Sz';
            int_cell{count+2, 4} = 'Sz';
            int_cell{count+2, 5} = Para.Model.JD;
            
            count = count + 3;
        end
    end
    
end

if strcmp(BCX, 'PBC')
    bnum = Ly * (Lx - 1);
    for y = 1:1:Ly
        
        int_cell{count, 1} = y;
        int_cell{count, 2} = bnum + y;
        int_cell{count, 3} = 'Sx';
        int_cell{count, 4} = 'Sx';
        int_cell{count, 5} = Para.Model.J;
        
        count = count + 1;
        
        int_cell{count, 1} = y;
        int_cell{count, 2} = bnum + y;
        int_cell{count, 3} = 'Sy';
        int_cell{count, 4} = 'Sy';
        int_cell{count, 5} = Para.Model.J;
        
        count = count + 1;
        
        int_cell{count, 1} = y;
        int_cell{count, 2} = bnum + y;
        int_cell{count, 3} = 'Sz';
        int_cell{count, 4} = 'Sz';
        int_cell{count, 5} = Para.Model.J;
        
        count = count + 1;
    end
    
    if mod(Lx, 2) == 0
        for y = 2:2:(Ly-1)
            % Jp: '/'
            int_cell{count, 1} = y;
            int_cell{count, 2} = bnum + y + 1;
            int_cell{count, 3} = 'Sx';
            int_cell{count, 4} = 'Sx';
            int_cell{count, 5} = Para.Model.JD;
            
            count = count + 1;
            
            int_cell{count, 1} = y;
            int_cell{count, 2} = bnum + y + 1;
            int_cell{count, 3} = 'Sy';
            int_cell{count, 4} = 'Sy';
            int_cell{count, 5} = Para.Model.JD;
            
            count = count + 1;
            
            int_cell{count, 1} = y;
            int_cell{count, 2} = bnum + y + 1;
            int_cell{count, 3} = 'Sz';
            int_cell{count, 4} = 'Sz';
            int_cell{count, 5} = Para.Model.JD;
            
            count = count + 1;
        end
    else
        for y = 2:2:Ly
            % Jp: '\'
            int_cell{count, 1} = y;
            int_cell{count, 2} = bnum + y - 1;
            int_cell{count, 3} = 'Sx';
            int_cell{count, 4} = 'Sx';
            int_cell{count, 5} = Para.Model.JD;
            
            count = count + 1;
            
            int_cell{count, 1} = y;
            int_cell{count, 2} = bnum + y - 1;
            int_cell{count, 3} = 'Sy';
            int_cell{count, 4} = 'Sy';
            int_cell{count, 5} = Para.Model.JD;
            
            count = count + 1;
            
            int_cell{count, 1} = y;
            int_cell{count, 2} = bnum + y - 1;
            int_cell{count, 3} = 'Sz';
            int_cell{count, 4} = 'Sz';
            int_cell{count, 5} = Para.Model.JD;
            
            count = count + 1;
        end
    end
end

if strcmp(BCY, 'PBC')
    for x = 1:1:Lx       
        
        int_cell{count, 1} = (x-1)*Ly+1;
        int_cell{count, 2} = x*Ly;
        int_cell{count, 3} = 'Sx';
        int_cell{count, 4} = 'Sx';
        int_cell{count, 5} = Para.Model.J;
        
        count = count + 1;
        
        int_cell{count, 1} = (x-1)*Ly+1;
        int_cell{count, 2} = x*Ly;
        int_cell{count, 3} = 'Sy';
        int_cell{count, 4} = 'Sy';
        int_cell{count, 5} = Para.Model.J;
        
        count = count + 1;
        
        int_cell{count, 1} = (x-1)*Ly+1;
        int_cell{count, 2} = x*Ly;
        int_cell{count, 3} = 'Sz';
        int_cell{count, 4} = 'Sz';
        int_cell{count, 5} = Para.Model.J;
        
        count = count + 1;
    end
    
    if mod(Ly, 2) == 0
        for x = 2:2:(Lx-1)
            
            int_cell{count, 1} = (x-1)*Ly+1;
            int_cell{count, 2} = x*Ly+Ly;
            int_cell{count, 3} = 'Sx';
            int_cell{count, 4} = 'Sx';
            int_cell{count, 5} = Para.Model.JD;
            
            count = count + 1;
            
            int_cell{count, 1} = (x-1)*Ly+1;
            int_cell{count, 2} = x*Ly+Ly;
            int_cell{count, 3} = 'Sy';
            int_cell{count, 4} = 'Sy';
            int_cell{count, 5} = Para.Model.JD;
            
            count = count + 1;
            
            int_cell{count, 1} = (x-1)*Ly+1;
            int_cell{count, 2} = x*Ly+Ly;
            int_cell{count, 3} = 'Sz';
            int_cell{count, 4} = 'Sz';
            int_cell{count, 5} = Para.Model.JD;
            
            count = count + 1;
        end
    else
        for x = 2:2:Lx
            
            int_cell{count, 1} = (x-1)*Ly+1;
            int_cell{count, 2} = x*Ly-Ly;
            int_cell{count, 3} = 'Sx';
            int_cell{count, 4} = 'Sx';
            int_cell{count, 5} = Para.Model.JD;
            
            count = count + 1;
            
            int_cell{count, 1} = (x-1)*Ly+1;
            int_cell{count, 2} = x*Ly-Ly;
            int_cell{count, 3} = 'Sy';
            int_cell{count, 4} = 'Sy';
            int_cell{count, 5} = Para.Model.JD;
            
            count = count + 1;
            
            int_cell{count, 1} = (x-1)*Ly+1;
            int_cell{count, 2} = x*Ly-Ly;
            int_cell{count, 3} = 'Sz';
            int_cell{count, 4} = 'Sz';
            int_cell{count, 5} = Para.Model.JD;
            
            count = count + 1;
        end
    end
end

Intr = struct('JmpOut', int_cell(:,1), 'JmpIn', int_cell(:,2), ...
              'JmpOut_type', int_cell(:,3), 'JmpIn_type', int_cell(:,4), 'CS', int_cell(:,5));

for i = 1:1:length(Intr)
    if Intr(i).JmpOut > Intr(i).JmpIn
        T = Intr(i).JmpIn;
        Intr(i).JmpIn = Intr(i).JmpOut;
        Intr(i).JmpOut = T;
    end
end
end