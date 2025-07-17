function [ PseudoMPO ] = Automata_SM( Para )
% function [ IntrMap ] = AutoMataInit( Para )
% Initialize the Hamiltonian by Automata picture.
% Yuan Gao@buaa 2021.04.06
% mail: 17231064@buaa.edu.cn
% keyboard;
IntrMap = eval([Para.IntrcMap_Name, '(Para)']);
% IntrMap = SortIntrMap( IntrMap );
L = Para.L;
d = Para.d;
BondInfo = cell(L, 1);

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

for i = 1:1:length(BondInfo)
    PseudoMPO(i).A = cell(BondInfo{i});
    if i ~= L
        PseudoMPO(i).A{1,1} = sparse(eye(d));
    end
    if i ~= 1
        PseudoMPO(i).A{BondInfo{i}(1),BondInfo{i}(2)} = sparse(eye(d));
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
        PseudoMPO(site_list(1)).A{1,end} = PseudoMPO(site_list(1)).A{1,end} + CS * sparse(op_list{1});
    else
        % two site --------------------
        pos = pos_list(site_list(1));
        PseudoMPO(site_list(1)).A{1, pos} = sparse(op_list{1});
        pos_list(site_list(1)) = pos_list(site_list(1)) + 1;
        
        if site_list(end) > site_list(1) + 1
            if length(site_list) == 2
                for j = (site_list(1)+1):1:(site_list(end)-1)
                    PseudoMPO(j).A{pos, pos_list(j)} = sparse(eye(d));
                    pos = pos_list(j);
                    pos_list(j) = pos_list(j) + 1;
                end
            else
                cnt = 2;
                for j = (site_list(1)+1):1:(site_list(end)-1)
                    if j == site_list(cnt)
                        PseudoMPO(j).A{pos, pos_list(j)} = sparse(op_list{cnt});
                        pos = pos_list(j);
                        pos_list(j) = pos_list(j) + 1;
                        cnt = cnt + 1;
                    else
                        PseudoMPO(j).A{pos, pos_list(j)} = sparse(eye(d));
                        pos = pos_list(j);
                        pos_list(j) = pos_list(j) + 1;
                    end
                end
            end
        end
        PseudoMPO(site_list(end)).A{pos,end} = CS * sparse(op_list{end});
    end
end


for it = 1:L
    szinfo = size(PseudoMPO(it).A);
    for iti = 1:szinfo(1)
        for itj = 1:szinfo(2)
            if isempty(PseudoMPO(it).A{iti, itj})
                PseudoMPO(it).A{iti, itj} = sparse(zeros([2,2]));
            end
        end
    end
end

end

