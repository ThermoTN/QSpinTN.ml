clear all
close all
L = 5; % system size
J1 = 1; % NN coupling
J2 = 0.5; % NNN coupling
h = 0.1; % magnetic field
d = 2; % dimension of local space

% operator
Sz = [-0.5, 0; 0, 0.5];
Sx = [0, 0.5; 0.5, 0];
Id = eye(2);

cnt = 1;
% NN interaction
for it = 1:(L-1)
    Intr(cnt).site = [it, it+1];
    Intr(cnt).operator = {Sz, Sz};
    Intr(cnt).J = J1;
    cnt = cnt + 1;
end

% NNN interaction
for it = 1:(L-2)
    Intr(cnt).site = [it, it+2];
    Intr(cnt).operator = {Sz, Sz};
    Intr(cnt).J = J2;
    cnt = cnt + 1;
end

% Magnetic field
for it = 1:L
    Intr(cnt).site = [it];
    Intr(cnt).operator = {Sx};
    Intr(cnt).J = -h;
    cnt = cnt + 1;
end


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
for i = 1:1:length(Intr)
    if length(Intr(i).site) > 1
        BondInfo{Intr(i).site(1)}(2) = BondInfo{Intr(i).site(1)}(2) + 1;
        BondInfo{Intr(i).site(end)}(1) = BondInfo{Intr(i).site(end)}(1) + 1;
        if Intr(i).site(end) > Intr(i).site(1)+1
            for j = (Intr(i).site(1)+1):1:(Intr(i).site(end)-1)
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
for i = 1:1:length(Intr)
    site_list = Intr(i).site;
    op_list = Intr(i).operator;
    CS = Intr(i).J;
    if length(site_list) == 1
        % on site ---------------------
        try
            PseudoMPO(site_list(1)).A{1,end} = PseudoMPO(site_list(1)).A{1,end} + CS * sparse(op_list{1});
        catch
            PseudoMPO(site_list(1)).A{1,end} = CS * sparse(op_list{1});
        end
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