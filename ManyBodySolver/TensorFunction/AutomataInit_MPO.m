function [ H ] = AutomataInit_MPO( IntrMap, Para, L )
% function [ IntrMap ] = AutoMataInit( Para )
% Initialize the Hamiltonian by Automata picture.
% Yuan Gao@buaa 2021.04.06 modified 2025.07.09
% mail: yuangao@buaa.edu.cn
d = Para.d;
BondInfo = cell(L, 1);
Op = GetLocalSpace(d);
for i = 1:1:L
    if i == 1
        BondInfo{i} = [1,2];
    elseif i == L
        BondInfo{i} = [2,1];
    else
        BondInfo{i} = [2,2];
    end
end

% =========== Get Bond Dim Info ===========
for i = 1:1:length(IntrMap)
    if IntrMap(i).site(1)>0
        BondInfo{IntrMap(i).site(1)}(2) = BondInfo{IntrMap(i).site(1)}(2) + 1;
        BondInfo{IntrMap(i).site(end)}(1) = BondInfo{IntrMap(i).site(end)}(1) + 1;
        if IntrMap(i).site(end) > IntrMap(i).site(1)+1
            for j = (IntrMap(i).site(1)+1):1:(IntrMap(i).site(end)-1)
                BondInfo{j} = BondInfo{j} + [1,1];
            end
        end
    end
end

% =========== Get MPO  ===========
H.A = cell(L,1);
for i = 1:1:length(BondInfo)
    H.A{i} = zeros([BondInfo{i},d,d]);
    %     PseudoMPO(i).A = cell(BondInfo{i});
    if i ~= L
        H.A{i}(1,1,:,:) = reshape(Op.Id, [1,1,d,d]);
        %         PseudoMPO(i).A{1,1} = 'Id';
    end
    if i ~= 1
        H.A{i}(BondInfo{i}(1),BondInfo{i}(2),:,:) = reshape(Op.Id, [1,1,d,d]);
        %         PseudoMPO(i).A{BondInfo{i}(1),BondInfo{i}(2)} = 'Id';
    end
end
pos_list = zeros(L, 1);
pos_list(:) = 2;
for i = 1:1:length(IntrMap)
    site_list = IntrMap(i).site;
    op_list = IntrMap(i).operator;
    CS = IntrMap(i).J;
    if length(site_list) == 1
        % on site ---------------------
        H.A{site_list(1)}(1,end,:,:) = CS * reshape(op_list{1}, [1,1,d,d]) + H.A{site_list(1)}(1,end,:,:);
    else
        % two site --------------------
        pos = pos_list(site_list(1));
        H.A{site_list(1)}(1,pos,:,:) = reshape(op_list{1}, [1,1,d,d]);
        pos_list(site_list(1)) = pos_list(site_list(1)) + 1;
        
        if site_list(end) > site_list(1) + 1
            if length(site_list) == 2
                for j = (site_list(1)+1):1:(site_list(end)-1)
                    H.A{j}(pos, pos_list(j),:,:) = reshape(Op.Id, [1,1,d,d]);
                    pos = pos_list(j);
                    pos_list(j) = pos_list(j) + 1;
                end
            else
                cnt = 2;
                for j = (site_list(1)+1):1:(site_list(end)-1)
                    if j == site_list(cnt)
                        H.A{j}(pos, pos_list(j),:,:) = reshape(op_list{cnt}, [1,1,d,d]);
                        pos = pos_list(j);
                        pos_list(j) = pos_list(j) + 1;
                        cnt = cnt + 1;
                    else
                        H.A{j}(pos, pos_list(j),:,:) = reshape(Op.Id, [1,1,d,d]);
                        pos = pos_list(j);
                        pos_list(j) = pos_list(j) + 1;
                    end
                end
            end
        end
        H.A{site_list(end)}(pos,end,:,:) = CS * reshape(op_list{end}, [1,1,d,d]);
    end
end
H.A{1} = permute(H.A{1}, [2,3,4,1]);
H.A{end} = permute(H.A{end}, [1,3,4,2]);
if length(size(H.A{1})) == 4
    H.A{1} = H.A{1}(:,:,:,1);
end

if length(size(H.A{end})) == 4
    H.A{end} = H.A{end}(:,:,:,end);
end
H.lgnorm = 0;
end

