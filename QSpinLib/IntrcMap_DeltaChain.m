function [ Intr ] = IntrcMap_DeltaChain( Para )
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
Lx = Para.Lx;
L = Lx;
h = Para.h;
int_num = L-1 + L/2-1 + L; % nn nnn on-site
int_num = int_num * 3;
int_cell = cell(int_num, 5);
count = 1;

% OBC part

% for i = 1:L
%     int_cell{count, 1} = 0;
%     int_cell{count, 2} = i;
%     int_cell{count, 3} = '~';
%     int_cell{count, 4} = 'Sx';
%     int_cell{count, 5} = h(1);
%
%     int_cell{count+1, 1} = 0;
%     int_cell{count+1, 2} = i;
%     int_cell{count+1, 3} = '~';
%     int_cell{count+1, 4} = 'Sy';
%     int_cell{count+1, 5} = h(2);
%
%     int_cell{count+2, 1} = 0;
%     int_cell{count+2, 2} = i;
%     int_cell{count+2, 3} = '~';
%     int_cell{count+2, 4} = 'Sz';
%     int_cell{count+2, 5} = h(3);
%
%     count = count + 3;
% end

for i = 1:L-1
    int_cell{count, 1} = i;
    int_cell{count, 2} = i+1;
    int_cell{count, 3} = 'Sx';
    int_cell{count, 4} = 'Sx';
    int_cell{count, 5} = Para.Model.J1xy;
    
    int_cell{count+1, 1} = i;
    int_cell{count+1, 2} = i+1;
    int_cell{count+1, 3} = 'Sy';
    int_cell{count+1, 4} = 'Sy';
    int_cell{count+1, 5} = Para.Model.J1xy;
    
    int_cell{count+2, 1} = i;
    int_cell{count+2, 2} = i+1;
    int_cell{count+2, 3} = 'Sz';
    int_cell{count+2, 4} = 'Sz';
    int_cell{count+2, 5} = Para.Model.J1z;
    
    count = count + 3;
end

for i = 1:2:L-3
    int_cell{count, 1} = i;
    int_cell{count, 2} = i+2;
    int_cell{count, 3} = 'Sx';
    int_cell{count, 4} = 'Sx';
    int_cell{count, 5} = Para.Model.J2xy;
    
    int_cell{count+1, 1} = i;
    int_cell{count+1, 2} = i+2;
    int_cell{count+1, 3} = 'Sy';
    int_cell{count+1, 4} = 'Sy';
    int_cell{count+1, 5} = Para.Model.J2xy;
    
    int_cell{count+2, 1} = i;
    int_cell{count+2, 2} = i+2;
    int_cell{count+2, 3} = 'Sz';
    int_cell{count+2, 4} = 'Sz';
    int_cell{count+2, 5} = Para.Model.J2z;
    
    count = count + 3;
end

% PBC
% int_cell{count, 1} = 1;
% int_cell{count, 2} = Para.L;
% int_cell{count, 3} = 'Sx';
% int_cell{count, 4} = 'Sx';
% int_cell{count, 5} = Para.Model.J1xy;
%
% int_cell{count+1, 1} = 1;
% int_cell{count+1, 2} = Para.L;
% int_cell{count+1, 3} = 'Sy';
% int_cell{count+1, 4} = 'Sy';
% int_cell{count+1, 5} = Para.Model.J1xy;
%
% int_cell{count+2, 1} = 1;
% int_cell{count+2, 2} = Para.L;
% int_cell{count+2, 3} = 'Sz';
% int_cell{count+2, 4} = 'Sz';
% int_cell{count+2, 5} = Para.Model.J1z;
%
% count = count + 3;
%
% int_cell{count, 1} = 1;
% int_cell{count, 2} = Para.L-1;
% int_cell{count, 3} = 'Sx';
% int_cell{count, 4} = 'Sx';
% int_cell{count, 5} = Para.Model.J2xy;
%
% int_cell{count+1, 1} = 1;
% int_cell{count+1, 2} = Para.L-1;
% int_cell{count+1, 3} = 'Sy';
% int_cell{count+1, 4} = 'Sy';
% int_cell{count+1, 5} = Para.Model.J2xy;
%
% int_cell{count+2, 1} = 1;
% int_cell{count+2, 2} = Para.L-1;
% int_cell{count+2, 3} = 'Sz';
% int_cell{count+2, 4} = 'Sz';
% int_cell{count+2, 5} = Para.Model.J2z;
%
% count = count + 3;


[Sx, Sy, Sz, Id] = SpinOp(Para.d);
for it = 1:length(int_cell)
    if ~isempty(int_cell{it, 1})
        Intr(it).site = sort([int_cell{it,1}, int_cell{it,2}]);
        Intr(it).operator = {eval(int_cell{it,3}), eval(int_cell{it,4})};
        Intr(it).J = int_cell{it,5};
    end
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