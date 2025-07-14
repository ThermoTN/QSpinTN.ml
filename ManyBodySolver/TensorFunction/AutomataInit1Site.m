function [ H, Ope, htilde ] = AutomataInit1Site( H, Op, Para )
% function [ H, Ope ] = AutoMataInit1Site( H, Op, Para )
% Initialize the Hamiltonian by Automata picture.
% Yuan Gao@buaa 2021.04.06
% mail: 17231064@buaa.edu.cn
d = Para.d;
L = Para.L;


if isfield(Para, 'g_fec')
    g_fec = Para.g_fec.g;
    g_dir = Para.g_fec.dir;
else
    g_fec = [1,1,1];
    g_dir{1} = [1,0,0];
    g_dir{2} = [0,1,0];
    g_dir{3} = [0,0,1];
    Para.Field.hdir = Para.Field.h / norm(Para.Field.h);
    Para.Field.h = norm(Para.Field.h);
end

Sop = cell(1,3);
htilde = zeros(1,3);
for i = 1:3
    gi_dir = g_dir{i};
    Sop{i} = gi_dir(1) * Op.Sx3 + gi_dir(2) * Op.Sy3 + gi_dir(3) * Op.Sz3;
    htilde(i) = g_fec(i) * sum(gi_dir .* Para.Field.hdir);
end

Sxtot.A = cell(L, 1);
Sxtot.lgnorm = 0;

Sytot.A = cell(L, 1);
Sytot.lgnorm = 0;

Sztot.A = cell(L, 1);
Sztot.lgnorm = 0;

Ope.Sxtot = GetAMMPO(Sxtot, Op.Id, Op.Sx);
Ope.Sytot = GetAMMPO(Sytot, Op.Id, Op.Sy);
Ope.Sztot = GetAMMPO(Sztot, Op.Id, Op.Sz);

h = Para.Field.h;
if norm(Para.Field.h) ~= 0
    Sm.A = cell(Para.L, 1);
    Sm.lgnorm = 0;
    SmOp = htilde(1) * Sop{1} + htilde(2) * Sop{2} + htilde(3) * Sop{3};
    hOp = -Para.Field.h * SmOp;
    Ope.Sm = GetAMMPO(Sm, Op.Id, SmOp);
    Ope.SmOp = SmOp;
    for i = 1:1:L
        if i == 1
            H.A{i}(end,:,:) = H.A{i}(end,:,:) + reshape(hOp, [1,d,d]);
        elseif i == L
            H.A{i}(1,:,:) = H.A{i}(1,:,:) + reshape(hOp, [1,d,d]);
        else
            H.A{i}(1, end,:,:) = H.A{i}(1, end,:,:) + reshape(hOp, [1,1,d,d]);
        end
    end
end


end

function [ MPO ] = GetAMMPO( MPO, Id, S )
d = length(Id(:,1));
for i = 1:1:length(MPO.A)
    if i == 1
        MPO.A{i} = zeros(1,2,d,d);
        MPO.A{i}(1,1,:,:) = reshape(Id, [1,1,d,d]);
        MPO.A{i}(1,2,:,:) = reshape(S, [1,1,d,d]);
    elseif i == length(MPO.A)
        MPO.A{i} = zeros(2,1,d,d);
        MPO.A{i}(2,1,:,:) = reshape(Id, [1,1,d,d]);
        MPO.A{i}(1,1,:,:) = reshape(S, [1,1,d,d]);
    else
        MPO.A{i} = zeros(2,2,d,d);
        MPO.A{i}(1,1,:,:) = reshape(Id, [1,1,d,d]);
        MPO.A{i}(2,2,:,:) = reshape(Id, [1,1,d,d]);
        MPO.A{i}(1,2,:,:) = reshape(S, [1,1,d,d]);
    end
end
MPO.A{1} = permute(MPO.A{1}, [2,3,4,1]);
MPO.A{end} = permute(MPO.A{end}, [1,3,4,2]);
end









