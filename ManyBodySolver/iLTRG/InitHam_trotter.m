function [ Hab, Hba ] = InitHam_trotter( t_ab, t_ba, Para )
% function [ Hab, Hba ] = InitHam_tortter( Para )

[Sx, Sy, Sz, Id] = SpinOp(Para.d);

Hab = 0;
Hba = 0;

for i = 1:1:length(t_ab(:,1))
    switch t_ab{i, 1}
        case 'Sx'
            Sl = Sx;
        case 'Sy'
            Sl = Sy;
        case 'Sz'
            Sl = Sz;
        case 'Sz2'
            Sl = Sz * Sz;
        case 'Id'
            Sl = Id;
    end
    
    switch t_ab{i, 2}
        case 'Sx'
            Sr = Sx;
        case 'Sy'
            Sr = Sy;
        case 'Sz'
            Sr = Sz;
        case 'Sz2'
            Sl = Sz * Sz;
        case 'Id'
            Sl = Id;
    end
    
    if norm(Sl-Sr) < 1e-10
        hab = t_ab{i,3} * real(kron(Sl, Sr));
    else
        hab = t_ab{i,3} * kron(Sl, Sr);
    end
    Hab = Hab + hab;
end

for i = 1:1:length(t_ba(:,1))
    switch t_ba{i, 1}
        case 'Sx'
            Sl = Sx;
        case 'Sy'
            Sl = Sy;
        case 'Sz'
            Sl = Sz;
        case 'Sz2'
            Sl = Sz * Sz;
        case 'Id'
            Sl = Id;
    end
    
    switch t_ba{i, 2}
        case 'Sx'
            Sr = Sx;
        case 'Sy'
            Sr = Sy;
        case 'Sz'
            Sr = Sz;
        case 'Sz2'
            Sl = Sz * Sz;
        case 'Id'
            Sl = Id;
    end
    
    if norm(Sl-Sr) < 1e-10
        hba = t_ba{i,3} * real(kron(Sl, Sr));
    else
        hba = t_ba{i,3} * kron(Sl, Sr);
    end
    Hba = Hba + hba;
end

% Hxa = kron(Sx, Id);
% Hxb = kron(Id, Sx);
% Hya = kron(Sy, Id);
% Hyb = kron(Id, Sy);
% Hza = kron(Sz, Id);
% Hzb = kron(Id, Sz);

if norm(Para.Field.h) ~= 0
    g_fec = Para.g_fec.g;
    g_dir = Para.g_fec.dir;
    Sop = cell(1,3);
    htilde = zeros(1,3);
    for i = 1:3
        gi_dir = g_dir{i};
        Sop{i} = gi_dir(1) * Sx + gi_dir(2) * Sy + gi_dir(3) * Sz;
        Sop{i} = kron(Sop{i}, Id) + kron(Id, Sop{i});
        htilde(i) = g_fec(i) * sum(gi_dir .* Para.Field.hdir);
        Sop{i} = Sop{i} * htilde(i);
    end
    
    
    
    Hf = - 1/2 * Para.Field.h * (Sop{1} + Sop{2} + Sop{3});
    Hab = Hab + Hf;
    Hba = Hba + Hf;
end
end

