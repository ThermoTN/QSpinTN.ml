function [T, En, EE, EnvDir, info] = DMRGMain(Para, H, initialState)
% Numeric 1/2 select fixed-site algorithms. "auto" switches on an energy rise.
% Optional initialState must be normalized and right canonical (center 1).
if nargin<3, initialState=[]; end
if ~isempty(initialState)
    assert(isstruct(initialState) && numel(initialState.A)==Para.L, ...
        'DMRG:InitialState','Initial MPS must contain Para.L tensors.');
    assert(all(cellfun(@(a)size(a,1)<=Para.MCrit && size(a,2)<=Para.MCrit, ...
        initialState.A(2:end-1))), 'DMRG:InitialState','Initial MPS exceeds MCrit.');
end
mode = "auto";
if isfield(Para,'DMRGUpdateSites'), mode = string(Para.DMRGUpdateSites); end
assert(isscalar(mode) && ismember(mode,["1","2","auto"]),'DMRGUpdateSites must be 1, 2 or "auto".');
validateattributes(Para.DMRGStepMax,{'numeric'},{'scalar','integer','positive','finite'});
env = DMRGEnvironmentStore(Para);
cleanup = onCleanup(@() env.finish());
EnvDir = env.Directory;
info = struct('Mode',mode,'TwoSiteSweeps',0,'OneSiteSweeps',0,'Converged',false, ...
    'StopReason',"sweep_limit",'SwitchSweep',[],'SwitchDeltaE',[], ...
    'SwitchEnergy',[],'OneSiteInitialEnergy',[], ...
    'TwoSiteSweepSeconds',[],'OneSiteSweepSeconds',[], ...
    'TwoSiteLocalIterations',[],'OneSiteLocalIterations',[]);
handoff = [];
if mode~="1"
    [T,En,EE,info.TwoSiteSweeps,reason,rise,stats] = ...
        twoSiteSweeps(Para,H,env,mode=="auto",initialState);
    info.TwoSiteSweepSeconds = stats.Seconds;
    info.TwoSiteLocalIterations = stats.Iterations;
    if reason~="switch"
        info.Converged = reason=="converged";
        info.StopReason = reason;
        fprintf('DMRG_STOP sites=2 sweeps=%d reason=%s\n',info.TwoSiteSweeps,reason);
        return
    end
    info.SwitchSweep = info.TwoSiteSweeps;
    info.SwitchDeltaE = rise;
    info.SwitchEnergy = real(En);
    handoff = struct('T',T,'Sweep',info.SwitchSweep,'Energy',info.SwitchEnergy);
    fprintf('DMRG_SWITCH 2->1 after sweep %d: E=%.16g, deltaE=%+.12g\n', ...
        info.SwitchSweep,info.SwitchEnergy,rise);
    env.reset();  % Discard fused two-site blocks before building single-site blocks.
elseif ~isempty(initialState)
    handoff = struct('T',initialState,'Sweep',0);
end
[T,En,EE,stage] = DMRGOneSite(Para,H,env,handoff);
info.OneSiteSweeps = stage.Sweeps;
info.OneSiteInitialEnergy = stage.InitialEnergy;
info.Converged = stage.Converged;
info.OneSiteSweepSeconds = stage.SweepSeconds;
info.OneSiteLocalIterations = stage.LocalIterations;
if stage.Converged, info.StopReason = "converged"; end
fprintf('DMRG_STOP sites=1 sweeps=%d total=%d reason=%s\n', ...
    stage.Sweeps,info.TwoSiteSweeps+stage.Sweeps,info.StopReason);
end

function [T,En,EE,i,reason,rise,stats] = twoSiteSweeps(Para,H,EnV,switchOnRise,initialState)
reason = "sweep_limit"; rise = [];
stats = struct('Seconds',[],'Iterations',[]);
sweepTol = GetDMRGSweepTol(Para);
if ~isempty(initialState)
    T = initialState;
else
switch Para.d
    case 2
        T = InitMPS(Para);
    otherwise
        T = InitMPS_d(Para);
end
end

Rslt.En = zeros(1, Para.DMRGStepMax+2);
Rslt.En(1) = inf;
% Rslt.En(2) = Observe( T, H, Para );


Init_DMRG_EnV(EnV,H,T);
Rslt.En(2) = GetEn_DMRG(EnV, T) * exp(H.lgnorm);
% keyboard;
EE = zeros([Para.L-1, 1]);
TE = zeros([Para.L-1, 1]);
maxEL2R = zeros([Para.L-1, 1]);
maxER2L = zeros([Para.L-1, 1]);
EEOld = 0;
DK = Para.DK;
for i = 1:1:Para.DMRGStepMax
    if i == 1 && isempty(initialState)
	Para.DK = 2;
    else
	Para.DK = DK;
    end
    ticSweep = tic; localIterations = 0;
    left = EnV.get(1,'L');
    for j = 1:1:(Para.L-2)
        [T, ~, maxEL2R(j), ~, k] = two_site_update_DMRG(left, EnV.get(j+1,'R'), T, j, '->-', Para);
        localIterations = localIterations+k;
        left = Update_DMRG_EnV(left, H, T, j, '->-');
        EnV.put(j+1,'L',left);
        % fprintf('  %d, ->-\n', j)
    end
    clear left
    right = EnV.get(Para.L,'R');
    for j = Para.L:-1:2
        [T, TE(j-1), maxER2L(j-1), ~, k] = two_site_update_DMRG(EnV.get(j-1,'L'), right, T, j, '-<-', Para);
        localIterations = localIterations+k;
        if j > 2
            right = Update_DMRG_EnV(right, H, T, j, '-<-');
            EnV.put(j-1,'R',right);
        end
        % fprintf('  %d, -<-\n', j)
    end
    clear right
    EE = max(maxER2L);
    % Rslt.En(i + 2) = Observe(T, H, Para);
    Rslt.En(i + 2) = GetEn_DMRG(EnV, T) * exp(H.lgnorm);
    [converged,DeltaS,deltaE] = CheckDMRGSweep( ...
        Rslt.En(i+2),Rslt.En(i+1),EE,EEOld,sweepTol);
    EEOld = EE;
    % keyboard;
    Rslt.MaxEE = EE;
    fprintf('%d, En: %.16f, DeltaS: %g\n     Trunc. Er.: %g\n     EE: %.16f\n', i, Rslt.En(i + 2), ...
        DeltaS, max(TE), EE);
    current = real(Rslt.En(i+2)); previous = real(Rslt.En(i+1));
    assert(isfinite(current),'DMRG:Energy','Nonfinite two-site energy at sweep %d.',i);
    stats.Seconds(i) = toc(ticSweep); stats.Iterations(i) = localIterations;
    fprintf('Elapsed time is %.6f seconds.\n',stats.Seconds(i));
    % Only compare completed sweeps. Give a resolved rise priority over stopping.
    roundoff = 16*eps(max([1,abs(current),abs(previous)]));
    if switchOnRise && i>=2 && deltaE>roundoff
        reason = "switch"; rise = deltaE; break
    end
    if i>=2 && converged
        reason = "converged"; break
    end
    if ~switchOnRise && i>=4 && deltaE>0
        reason = "energy_rise"; break
    end
end
En = Rslt.En(i + 2);
end

% ======================================== %
% Sub Function
% ======================================== %
function Init_DMRG_EnV(EnV,H,T)
EnV.put(1,'L',H.A{1});
EnV.put(numel(H.A),'R',H.A{end});
Vr = H.A{end};

for i = (length(H.A)-1):-1:2
    if i == (length(H.A)-1)
        Vr = contract(T.A{i+1}, 2, Vr, 3);
        Vr = contract(conj(T.A{i+1}), 2, Vr, 3);
        Vr = contract(Vr, 3, H.A{i},2);
    else
        Vr = contract(T.A{i+1}, [2,3], Vr, [2,5]);
        Vr = contract(conj(T.A{i+1}), [2,3], Vr, [2,4]);
        Vr = contract(Vr, 3, H.A{i}, 2);
    end
    EnV.put(i,'R',Vr);
end
end

function [T, TE, EE, maxE, iterations] = two_site_update_DMRG(Vl, Vr, T, pos, dir, Para)

switch dir
    case '->-'
        if pos == 1
            TT = contract(T.A{pos}, 1, T.A{pos+1}, 1);
        else
            TT = contract(T.A{pos}, 2, T.A{pos+1}, 1);
        end
    case '-<-'
        if pos == 2
            TT = contract(T.A{pos-1}, 1, T.A{pos}, 1);
        else
            TT = contract(T.A{pos-1}, 2, T.A{pos}, 1);
        end
end
% Both modes use the same three-term recurrence and local stopping rule.
[TT,~,iterations] = DMRGLanczos(@(x) twoSiteAction(x,Vl,Vr),TT,Para);

Tsize = size(TT);
if pos == 1 || (pos == 2 && strcmp(dir, '-<-'))
    TT = reshape(TT, [Tsize(1), prod(Tsize(2:1:3))]);
elseif (pos == (Para.L-1) && strcmp(dir, '->-')) || pos == Para.L
    TT = reshape(TT, [prod(Tsize(1:1:2)), prod(Tsize(3))]);
else
    TT = reshape(TT, [prod(Tsize(1:1:2)), prod(Tsize(3:1:4))]);
end

[U, S, V, TE, EE] = svdT(TT, 'Nkeep', Para.MCrit, 'epsilon', 1e-8);
Usize = size(U);
Vsize = size(V);
if pos == 1 || (pos == 2 && strcmp(dir, '-<-'))
    U = permute(U, [2,1]);
    V = reshape(V, [Vsize(1), Tsize(2), Tsize(3)]);
elseif (pos == (Para.L-1) && strcmp(dir, '->-')) || pos == Para.L
    U = reshape(U, [Tsize(1), Tsize(2), Usize(2)]);
    U = permute(U, [1,3,2]);
else
    U = reshape(U, [Tsize(1), Tsize(2), Usize(2)]);
    U = permute(U, [1,3,2]);
    V = reshape(V, [Vsize(1), Tsize(3), Tsize(4)]);
end

Ns = norm(S);
S = diag(S)/Ns;
maxE = max(diag(S));
switch dir
    case '->-'
        T.A{pos} = U;
        T.A{pos+1} = contract(S, 2, V, 1);
    case '-<-'
        if pos-1 == 1
            T.A{pos-1} = contract(S, 1, U, 1);
        else
            T.A{pos-1} = contract(U, 2, S, 1, [1,3,2]);
        end
        T.A{pos} = V;
end
%keyboard;
end

function Y = twoSiteAction(X,Vl,Vr)
if ndims(Vl)==3
    Y = contract(X,1,Vl,3);
else
    Y = contract(X,[1,2],Vl,[2,5]);
    if ndims(Y)==5, Y=permute(Y,[1,2,4,3,5]);
    else, Y=permute(Y,[1,3,2,4]); end
end
if ndims(Vr)==3
    Y = contract(Y,[1,2],Vr,[3,1]);
else
    Y = contract(Y,[1,2,3],Vr,[2,5,3]);
end
end

function block = Update_DMRG_EnV(block, H, T, pos, dir)
switch dir
    case '->-'
        if pos == 1
            block = contract(T.A{pos}, 2, block, 3);
            block = contract(conj(T.A{pos}), 2, block, 3);
            block = contract(block, 3, H.A{pos+1}, 1);
        else
            block = contract(T.A{pos}, [1,3], block, [2,5]);
            block = contract(conj(T.A{pos}), [1,3], block, [2,4]);
            block = contract(block, 3, H.A{pos+1}, 1);
        end
    case '-<-'
        if pos == length(H.A)
            block = contract(T.A{pos}, 2, block, 3);
            block = contract(conj(T.A{pos}), 2, block, 3);
            block = contract(block, 3, H.A{pos-1}, 2);
        else
            block = contract(T.A{pos}, [2,3], block, [2,5]);
            block = contract(conj(T.A{pos}), [2,3], block, [2,4]);
            block = contract(block, 3, H.A{pos-1}, 2);
        end
end
end
