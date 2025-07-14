function [ Intr ] = IntrcMap_XXZtest( Para )
% Usage: [ Intr ] = IntrcMap_XXZtest( Para ) 
% 1D spin chain XXZ model
% Model parameter:
%           Para.Model.Jxy      Jxy     Nearest neighbor xy term
%           Para.Model.Jz       Jz      Nearest neighbor z term
%
% Lattice geometry:
%           Para.Geo.L          Lattice size
%           Para.Geo.BC         Boundary condition
%                                   'PBC', periodic boundary condition (add S_1 S_L term)
%                                   'OBC', open boundary condition
%
% Hamiltonian (OBC):
%       H = sum_i^L-1 Jxy (Sx_i Sx_i+1 + Sy_i Sy_i+1) + Jz (Sz_i Sz_i+1)

L = Para.Geo.L;
if isinf(L)
    error('Illegal lattice sites!')
end
int_num = 3 * (L-1);
BC = Para.Geo.BC;

if strcmp(BC, 'PBC')
    int_num = int_num + 3;
end
int_cell = cell(int_num, 5);
count = 1;

% OBC part

for x = 1:1:(L-1)
    
    int_cell{count, 1} = x;
    int_cell{count, 2} = x + 1;
    int_cell{count, 3} = 'Sx';
    int_cell{count, 4} = 'Sx';
    int_cell{count, 5} = Para.Model.Jxy;
    count = count + 1;
    
    int_cell{count, 1} = x;
    int_cell{count, 2} = x + 1;
    int_cell{count, 3} = 'Sy';
    int_cell{count, 4} = 'Sy';
    int_cell{count, 5} = Para.Model.Jxy;
    count = count + 1;
    
    int_cell{count, 1} = x;
    int_cell{count, 2} = x + 1;
    int_cell{count, 3} = 'Sz';
    int_cell{count, 4} = 'Sz';
    int_cell{count, 5} = Para.Model.Jz;
    
    count = count + 1;
end

if strcmp(BC, 'PBC')
    int_cell{count, 1} = 1;
    int_cell{count, 2} = L;
    int_cell{count, 3} = 'Sx';
    int_cell{count, 4} = 'Sx';
    int_cell{count, 5} = Para.Model.Jxy;
    count = count + 1;
    
    int_cell{count, 1} = 1;
    int_cell{count, 2} = L;
    int_cell{count, 3} = 'Sy';
    int_cell{count, 4} = 'Sy';
    int_cell{count, 5} = Para.Model.Jxy;
    count = count + 1;
    
    int_cell{count, 1} = 1;
    int_cell{count, 2} = L;
    int_cell{count, 3} = 'Sz';
    int_cell{count, 4} = 'Sz';
    int_cell{count, 5} = Para.Model.Jz;
    
    count = count + 1;
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