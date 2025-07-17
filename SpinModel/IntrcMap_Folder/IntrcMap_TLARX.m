function [ Intr ] = IntrcMap_TLARX( Para )
% 2022.4.14 yGao & ezLv
% Triangular lattice
% XXZ model
% Parameter:
%           J1xy        Nearest neighbor SxSx+SySy term
%           J1z         Nearest neighbor SzSz term
%           J2xy        Next-nearest neighbor SxSx+SySy term
%           J2z         Next-nearest neighbor SzSz term
%
% Hamiltonian:
%   H = \sum_<i,j> J1xy (Sx_i Sx_j + Sy_i Sy_j) + J1z Sz_i Sz_j
%       + \sum_<<i,j>> J1xy (Sx_i Sx_j + Sy_i Sy_j) + J1z Sz_i Sz_j
%       - h\sum_i Sh_i
% keyboard;
[ Bond_info, ~, ~ ] = Bond_Info_TLARX( Para );
% keyboard;

int_cell = cell(1, 5);
count = 1;

% OBC part

for i = 1:1:length(Bond_info)
    [int_cell, count] = add_intr(Para, int_cell, Bond_info(i, :), count);
end

[Sx, Sy, Sz, Id] = SpinOp(Para.d);
for it = 1:length(int_cell)
    [Intr(it).site, ord] = sort([int_cell{it,1}, int_cell{it,2}]);
    if ord(1) == 2
        Intr(it).operator = {eval(int_cell{it,4}), eval(int_cell{it,3})};
    else
        Intr(it).operator = {eval(int_cell{it,3}), eval(int_cell{it,4})};
    end
    Intr(it).J = int_cell{it,5};
    
end

% Intr = struct('JmpOut', int_cell(:,1), 'JmpIn', int_cell(:,2), ...
%     'JmpOut_type', int_cell(:,3), 'JmpIn_type', int_cell(:,4), 'CS', int_cell(:,5));
% 
% for i = 1:1:length(Intr)
%     if Intr(i).JmpOut > Intr(i).JmpIn
%         T = Intr(i).JmpIn;
%         Intr(i).JmpIn = Intr(i).JmpOut;
%         Intr(i).JmpOut = T;
%     end
% end
end

function [int_cell, n] = add_intr(Para, int_cell, Bond_info, n)

JmpOut = Bond_info{1};
JmpIn = Bond_info{2};
Bond_type = Bond_info{3};
switch Bond_type(1:1:3)
    case 'nn-'
        
        % xxz term
        int_cell{n, 1} = JmpOut;
        int_cell{n, 2} = JmpIn;
        int_cell{n, 3} = 'Sx';
        int_cell{n, 4} = 'Sx';
        int_cell{n, 5} = Para.Model.J1xy;
        
        int_cell{n+1, 1} = JmpOut;
        int_cell{n+1, 2} = JmpIn;
        int_cell{n+1, 3} = 'Sy';
        int_cell{n+1, 4} = 'Sy';
        int_cell{n+1, 5} = Para.Model.J1xy;
        
        int_cell{n+2, 1} = JmpOut;
        int_cell{n+2, 2} = JmpIn;
        int_cell{n+2, 3} = 'Sz';
        int_cell{n+2, 4} = 'Sz';
        int_cell{n+2, 5} = Para.Model.J1z;
        
        n = n + 3;
        switch Bond_type(end)
            case '1'
                
                % JPD term
                int_cell{n, 1} = JmpOut;
                int_cell{n, 2} = JmpIn;
                int_cell{n, 3} = 'Sx';
                int_cell{n, 4} = 'Sx';
                int_cell{n, 5} = Para.Model.JPD * 2;
                
                int_cell{n+1, 1} = JmpOut;
                int_cell{n+1, 2} = JmpIn;
                int_cell{n+1, 3} = 'Sy';
                int_cell{n+1, 4} = 'Sy';
                int_cell{n+1, 5} = -Para.Model.JPD * 2;
                
                n = n + 2;
                
                % JGamma term
                int_cell{n, 1} = JmpOut;
                int_cell{n, 2} = JmpIn;
                int_cell{n, 3} = 'Sz';
                int_cell{n, 4} = 'Sy';
                int_cell{n, 5} = Para.Model.JGamma;
                
                int_cell{n+1, 1} = JmpOut;
                int_cell{n+1, 2} = JmpIn;
                int_cell{n+1, 3} = 'Sy';
                int_cell{n+1, 4} = 'Sz';
                int_cell{n+1, 5} = Para.Model.JGamma;
                
                n = n + 2;
            case '2'
                
                % JPD term
                int_cell{n, 1} = JmpOut;
                int_cell{n, 2} = JmpIn;
                int_cell{n, 3} = 'Sx';
                int_cell{n, 4} = 'Sx';
                int_cell{n, 5} = -Para.Model.JPD;
                
                int_cell{n+1, 1} = JmpOut;
                int_cell{n+1, 2} = JmpIn;
                int_cell{n+1, 3} = 'Sx';
                int_cell{n+1, 4} = 'Sy';
                int_cell{n+1, 5} = -Para.Model.JPD*sqrt(3);
                
                int_cell{n+2, 1} = JmpOut;
                int_cell{n+2, 2} = JmpIn;
                int_cell{n+2, 3} = 'Sy';
                int_cell{n+2, 4} = 'Sx';
                int_cell{n+2, 5} = -Para.Model.JPD*sqrt(3);
                
                int_cell{n+3, 1} = JmpOut;
                int_cell{n+3, 2} = JmpIn;
                int_cell{n+3, 3} = 'Sy';
                int_cell{n+3, 4} = 'Sy';
                int_cell{n+3, 5} = Para.Model.JPD;
                n = n + 4;
                
                % JGamma term
                int_cell{n, 1} = JmpOut;
                int_cell{n, 2} = JmpIn;
                int_cell{n, 3} = 'Sz';
                int_cell{n, 4} = 'Sx';
                int_cell{n, 5} = -Para.Model.JGamma*sqrt(3)/2;
                
                int_cell{n+1, 1} = JmpOut;
                int_cell{n+1, 2} = JmpIn;
                int_cell{n+1, 3} = 'Sx';
                int_cell{n+1, 4} = 'Sz';
                int_cell{n+1, 5} = -Para.Model.JGamma*sqrt(3)/2;
                
                int_cell{n+2, 1} = JmpOut;
                int_cell{n+2, 2} = JmpIn;
                int_cell{n+2, 3} = 'Sz';
                int_cell{n+2, 4} = 'Sy';
                int_cell{n+2, 5} = -Para.Model.JGamma/2;
                
                int_cell{n+3, 1} = JmpOut;
                int_cell{n+3, 2} = JmpIn;
                int_cell{n+3, 3} = 'Sy';
                int_cell{n+3, 4} = 'Sz';
                int_cell{n+3, 5} = -Para.Model.JGamma/2;
                
                n = n + 4;
            case '3'
                
                % JPD term
                int_cell{n, 1} = JmpOut;
                int_cell{n, 2} = JmpIn;
                int_cell{n, 3} = 'Sx';
                int_cell{n, 4} = 'Sx';
                int_cell{n, 5} = -Para.Model.JPD;
                
                int_cell{n+1, 1} = JmpOut;
                int_cell{n+1, 2} = JmpIn;
                int_cell{n+1, 3} = 'Sx';
                int_cell{n+1, 4} = 'Sy';
                int_cell{n+1, 5} = Para.Model.JPD*sqrt(3);
                
                int_cell{n+2, 1} = JmpOut;
                int_cell{n+2, 2} = JmpIn;
                int_cell{n+2, 3} = 'Sy';
                int_cell{n+2, 4} = 'Sx';
                int_cell{n+2, 5} = Para.Model.JPD*sqrt(3);
                
                int_cell{n+3, 1} = JmpOut;
                int_cell{n+3, 2} = JmpIn;
                int_cell{n+3, 3} = 'Sy';
                int_cell{n+3, 4} = 'Sy';
                int_cell{n+3, 5} = Para.Model.JPD;
                n = n + 4;
                
                % JGamma term
                int_cell{n, 1} = JmpOut;
                int_cell{n, 2} = JmpIn;
                int_cell{n, 3} = 'Sz';
                int_cell{n, 4} = 'Sx';
                int_cell{n, 5} = Para.Model.JGamma*sqrt(3)/2;
                
                int_cell{n+1, 1} = JmpOut;
                int_cell{n+1, 2} = JmpIn;
                int_cell{n+1, 3} = 'Sx';
                int_cell{n+1, 4} = 'Sz';
                int_cell{n+1, 5} = Para.Model.JGamma*sqrt(3)/2;
                
                int_cell{n+2, 1} = JmpOut;
                int_cell{n+2, 2} = JmpIn;
                int_cell{n+2, 3} = 'Sz';
                int_cell{n+2, 4} = 'Sy';
                int_cell{n+2, 5} = -Para.Model.JGamma/2;
                
                int_cell{n+3, 1} = JmpOut;
                int_cell{n+3, 2} = JmpIn;
                int_cell{n+3, 3} = 'Sy';
                int_cell{n+3, 4} = 'Sz';
                int_cell{n+3, 5} = -Para.Model.JGamma/2;
                
                n = n + 4;
        end
    case 'nnn'
        
        % xxz term
        int_cell{n, 1} = JmpOut;
        int_cell{n, 2} = JmpIn;
        int_cell{n, 3} = 'Sx';
        int_cell{n, 4} = 'Sx';
        int_cell{n, 5} = Para.Model.J2xy;
        
        int_cell{n+1, 1} = JmpOut;
        int_cell{n+1, 2} = JmpIn;
        int_cell{n+1, 3} = 'Sy';
        int_cell{n+1, 4} = 'Sy';
        int_cell{n+1, 5} = Para.Model.J2xy;
        
        int_cell{n+2, 1} = JmpOut;
        int_cell{n+2, 2} = JmpIn;
        int_cell{n+2, 3} = 'Sz';
        int_cell{n+2, 4} = 'Sz';
        int_cell{n+2, 5} = Para.Model.J2z;
        
        n = n + 3;
end
end



