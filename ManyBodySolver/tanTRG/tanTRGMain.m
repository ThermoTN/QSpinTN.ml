function Rslt = tanTRGMain(Para,H,rho,Op)
% Evolve a right-canonical, Frobenius-normalized thermal MPO from beta=2*tau0.
mode="auto";
if isfield(Para,'TanTRGUpdateSites'), mode=string(Para.TanTRGUpdateSites); end
assert(isscalar(mode) && ismember(mode,["1","2","cbe","auto"]), ...
    'TanTRG:Mode','TanTRGUpdateSites must be 1, 2, "cbe" or "auto".');
env=InitTSRGEnV(rho,H,Para); cleanup=onCleanup(@()env.finish());
Rslt=struct('beta',2*Para.tau0,'LnZ',2*rho.lgnorm/Para.L,'M',[], ...
    'En',real(GetEn(env.get(2,'R'),H,rho))*exp(H.lgnorm)/Para.L, ...
    'TrunErr',0,'EE',0);
if ~isempty(env.Directory), Rslt.EnvDir=env.Directory; end
Rslt=caltdvpThDQ(Para,rho,Op,Rslt);
Rslt.TanTRGInfo=struct('Mode',mode,'Sweeps',[]);
active=mode; if mode=="auto", active="2"; end
while Rslt.beta(end)<Para.beta_max
    Para.tau=GetNowtau(Para,Rslt.beta(end));
    validateattributes(Para.tau,{'numeric'},{'scalar','positive','finite'});
    Para.TanTRGUpdateSites=active;
    [rho,info]=TanTRGSweep(rho,H,Para,env);
    Rslt.beta(end+1)=Rslt.beta(end)+4*Para.tau;
    Rslt.LnZ(end+1)=2*rho.lgnorm/Para.L;
    Rslt.En(end+1)=real(GetEn(env.get(2,'R'),H,rho))*exp(H.lgnorm)/Para.L;
    Rslt.TrunErr(end+1)=info.TruncationError; Rslt.EE(end+1)=info.Entropy;
    info.Beta=Rslt.beta(end); Rslt.TanTRGInfo.Sweeps=[Rslt.TanTRGInfo.Sweeps,info];
    Rslt=caltdvpThDQ(Para,rho,Op,Rslt);
    fprintf(Para.fileID,'tanTRG sites=%s beta=%.16g En=%.16g D=%d time=%.6f s\n', ...
        active,Rslt.beta(end),Rslt.En(end),info.MaxBond,info.SweepSeconds);
    % Optional hooks are detected; errors in installed hooks must remain visible.
    if exist('GettanTRGCorr','file')==2, Rslt=GettanTRGCorr(rho,Rslt,Para,Op); end
    if exist('savetanTRGRslt','file')==2, savetanTRGRslt(Para,Rslt); end
    if mode=="auto" && ~Para.Full2 && Rslt.beta(end)>Para.beta_switch, active="1"; end
end
Rslt.rho=rho;
end
