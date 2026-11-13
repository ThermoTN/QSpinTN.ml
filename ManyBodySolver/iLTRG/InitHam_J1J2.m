function [Hab,Hba] = InitHam_J1J2(t_ab,t_ba,t_nnn,Para)
%INITHAM_J1J2 Spin-1/2 NN/NNN chain represented by two-spin blocks.
% Tables contain dimensionless coefficients and Sx/Sy/Sz operator names.
% Within block k=(A_k,B_k), t_ab is A_k -> B_k. t_ba connects
% B_k -> A_{k+1}; t_nnn connects A_k -> A_{k+1} and B_k -> B_{k+1}.
% Each block on-site Hamiltonian is shared equally by its two block bonds.
% Para.Field.h is the dimensionless Zeeman energy, with g already included.
% Copyright (c) 2025 ThermoTN. MIT License; see repository LICENSE.
validateattributes(Para.Field.h,{'numeric'},{'vector','numel',3,'real','finite'});
[Sx,Sy,Sz]=SpinOp(2);
spin=struct('Sx',Sx,'Sy',Sy,'Sz',Sz);
names={'Sx','Sy','Sz'}; left=struct(); right=struct();
for k=1:3
    name=names{k};
    left.(name)=kron(spin.(name),eye(2));
    right.(name)=kron(eye(2),spin.(name));
end
validate_table(t_ab,names); validate_table(t_ba,names); validate_table(t_nnn,names);
H0=zeros(4);
for k=1:size(t_ab,1)
    H0=H0+t_ab{k,3}*kron(spin.(t_ab{k,1}),spin.(t_ab{k,2}));
end
for k=1:3
    H0=H0-Para.Field.h(k)*(left.(names{k})+right.(names{k}));
end
Hab=(kron(H0,eye(4))+kron(eye(4),H0))/2;
for k=1:size(t_ba,1)
    Hab=Hab+t_ba{k,3}*kron(right.(t_ba{k,1}),left.(t_ba{k,2}));
end
for k=1:size(t_nnn,1)
    a=t_nnn{k,1}; b=t_nnn{k,2};
    Hab=Hab+t_nnn{k,3}*(kron(left.(a),left.(b))+kron(right.(a),right.(b)));
end
% Both tensor bonds repeat the same physical block orientation. Do not
% transpose or reverse directed DM terms on the second tensor bond.
Hba=Hab;
end

function validate_table(terms,names)
if ~iscell(terms) || size(terms,2)~=3
    error('qspintn:iltrg:Terms','Expected an operator/operator/coefficient table.');
end
for k=1:size(terms,1)
    if ~all(ismember(terms(k,1:2),names))
        error('qspintn:iltrg:Terms','Blocked spin-1/2 terms require Sx/Sy/Sz operators.');
    end
    validateattributes(terms{k,3},{'numeric'},{'scalar','real','finite'});
end
end
