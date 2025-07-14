clear all
format long
Rslt_data = [[0.110434143715077,0.00249778871224733,0.00249778871040138,0.00169861187929347,0.000389573624616153,0.00169861187641722,0.000316492751102309,0.000316492746969032,0.000416596157667703;0,0.110434143733950,0.00249778871208699,0.000389573624762706,0.00169861188090881,0.00169861187374651,0.000416596154271079,0.000316492745208230,0.000316492749356271;0,0,0.110434143778093,0.00169861187499778,0.00169861187933537,0.000389573612261552,0.000316492749705372,0.000416596152927552,0.000316492749795698;0,0,0,0.0554635046028409,0.00105023137365643,0.00105023137545706,0.000389573612618010,0.00169861187758101,0.00169861186851971;0,0,0,0,0.0554635045655813,0.00105023137844185,0.00169861187033708,0.00169861187071732,0.000389573617950223;0,0,0,0,0,0.0554635046350778,0.00169861187026119,0.000389573618617584,0.00169861186143272;0,0,0,0,0,0,0.110434143743894,0.00249778871781967,0.00249778871518136;0,0,0,0,0,0,0,0.110434143757569,0.00249778871523519;0,0,0,0,0,0,0,0,0.110434143714448]];
% addpath CommonQSfunc/
Sx = [0,1/2
    1/2,0];

Sy = [0,-1i/2
    1i/2,0];

Sz = [1/2,0
    0,-1/2];

Id = [1,0
    0,1];

% H = kron(Sz, Sz) %+ kron(Sy, Sy)
% keyboard;

L = 8;
Sztot = 0;
for i = 1:L
    Sx_list{i} = getOP(Sx, Id, i, L);
    Sy_list{i} = getOP(Sy, Id, i, L);
    Sz_list{i} = getOP(Sz, Id, i, L);
    Sp_list{i} = (Sx_list{i} + 1i * Sy_list{i});
    Sm_list{i} = Sp_list{i}';
    Sztot = Sztot + Sz_list{i};
end

H = 0;
       
% 
%   1-  4-  7
%    \ / \ / \
%     2-  5-  8
%      \ / \ / \
%       3-  6-  9
%        \ / \ / \
       
% Intr = {[1,2],[2,3],[4,5],[5,6],[7,8],[8,9],[1,4],[2,5],[3,6],[4,7],[5,8],[6,9],[1,3],[4,6],[7,9],[1,7],[2,8],[3,9]};
% for it = 1:length(Intr)
%     H = H + addterm(Sx_list, Sy_list, Sz_list, Intr{it}, -[1,1,1]);
% end
J = [1,1,1];
for it = 1:(L-1)
    H = H + addterm(Sx_list, Sy_list, Sz_list, [it, it + 1], J);
end

H0 = H;
for it = 1:L
    H = H - sqrt(3) * Sz_list{it};
end
beta = 16;
Z = trace(expm(-beta * H));
En = 1/Z * trace(expm(-beta * H) * H)/8;
load BenchMark_SpinChain.mat
loc = 8;
(trace(expm(-beta/2 * H) * (H * Sm_list{loc} - Sm_list{loc} * H) * expm(-beta/2 * H) * Sp_list{loc})/Z ...
    - (Rslt.ObsJ{end}(loc) + Rslt.ObsB{end}(loc))) ./ (Rslt.ObsJ{end}(loc) + Rslt.ObsB{end}(loc))
% 


% loc = 3;
% beta = 50;
% trace(expm(-beta/2 * H) * Sx_list{loc} * expm(-beta/2 * H) * Sx_list{loc})/trace(expm(-beta * H))
% Sztot = 0;
% for it = 1:L
%     Sztot = Sztot + Sz_list{it};
% end
% trace(expm(-beta * H) * Sztot * Sztot)/trace(expm(-beta * H))
% trace(expm(-beta/2 * H) * Sz_list{1} * expm(-beta/2 * H) * Sz_list{3})/trace(expm(-beta * H)) * 81
% Intr = {[01,03], [01,04], [04,06], [06,09], [07,09]};
% for it = 1:length(Intr)
%     H = H + addterm(Sx_list, Sy_list, Sz_list, Intr{it}, [0.5,0.5,1]);
% end
% h = 0.5;
% H = H - h * Sztot;
% h = 0;
% Para.Lx = 3;
% Para.Ly = 3;
% Sztot = 0;
% for it = 1:3
%     H = H - h * Sz_list{it};
%     Sztot = Sztot + Sz_list{it};
% end
% [V, D] = eig(H);
% GS = V(:,1);
% GS' * Sztot * GS/9
% beta = 16;
% Sm2 = Sm_list{1};
% Sm1 = Sm_list{1};
% Z = trace(expm(-beta * H));
% trace(expm(-beta/2 * H) * (H * Sm1 - Sm1 * H) * expm(-beta/2 * H) * Sm2')/Z



function Op = getOP(Op, Id, loc, L)
for i = 1:loc-1
    Op = kron(Id, Op);
end
for i = loc+1:L
    Op = kron(Op, Id);
end
end

function H = addterm(Sx_list, Sy_list, Sz_list, loc, J)
H = J(1) * Sx_list{loc(1)} * Sx_list{loc(2)} + J(2) * Sy_list{loc(1)} * Sy_list{loc(2)} + J(3) * (Sz_list{loc(1)} * Sz_list{loc(2)});
end