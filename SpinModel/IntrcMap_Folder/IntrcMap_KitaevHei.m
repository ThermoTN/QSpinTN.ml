function [Intr] = IntrcMap_KitaevHei(Para)

[ Bond_info, ~, ~ ] = Bond_Info_Kitaev( Para );
% keyboard;
% keyboard;
int_cell = cell(1, 5);
count = 1;

% OBC part

for i = 1:1:length(Bond_info)
    [int_cell, count] = add_intr(Para, int_cell, Bond_info(i, :), count);
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

% Intr = struct('JmpOut', int_cell(:,1), 'JmpIn', int_cell(:,2), ...
%               'JmpOut_type', int_cell(:,3), 'JmpIn_type', int_cell(:,4), 'CS', int_cell(:,5));
% 
% for i = 1:1:length(Intr)
%     if Intr(i).JmpOut > Intr(i).JmpIn
%         T = Intr(i).JmpIn;
%         Intr(i).JmpIn = Intr(i).JmpOut;
%         Intr(i).JmpOut = T;
%     end
% end
end

function [int_cell, count] = add_intr(Para, int_cell, Bond_info, count)
int_cell{count, 1} = Bond_info{1};
int_cell{count, 2} = Bond_info{2};
switch Bond_info{3}
    case 'x-type'
        int_cell{count, 3} = 'Sx';
        int_cell{count, 4} = 'Sx';
    case 'y-type'
        int_cell{count, 3} = 'Sy';
        int_cell{count, 4} = 'Sy';
    case 'z-type'
        int_cell{count, 3} = 'Sz';
        int_cell{count, 4} = 'Sz';
end
int_cell{count, 5} = Para.Model.K;

count = count + 1;
end