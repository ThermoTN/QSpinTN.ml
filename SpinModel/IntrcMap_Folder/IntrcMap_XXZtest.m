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

[Sx, Sy, Sz, Id] = SpinOp(Para.d);
Sy = imag(Sy);
for it = 1:length(int_cell)
    Intr(it).site = sort([int_cell{it,1}, int_cell{it,2}]);
    Intr(it).operator = {eval(int_cell{it,3}), eval(int_cell{it,4})};
    if strcmp(int_cell{it,3}, 'Sy')
        Intr(it).J = -int_cell{it,5};
    else
        Intr(it).J = int_cell{it,5};
    end
end

end