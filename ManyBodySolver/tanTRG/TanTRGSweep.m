function [rho,info] = TanTRGSweep(rho,H,Para,env)
% Exactly one forward/backward sweep, delta beta=4*Para.tau.
% Optional caller-owned env avoids rebuilding environments between sweeps.
mode=string(Para.TanTRGUpdateSites);
assert(isscalar(mode) && ismember(mode,["1","2","cbe"]), ...
    'TanTRG:Mode','TanTRGSweep requires TanTRGUpdateSites=1, 2 or "cbe".');
validateattributes(Para.tau,{'numeric'},{'scalar','positive','finite'});
if nargin<4
    env=InitTSRGEnV(rho,H,Para); cleanup=onCleanup(@()env.finish());
end
timer=tic;
if mode=="2"
    [rho,info]=TwoSiteTanTRGSweep(rho,H,Para,env);
else
    [rho,info]=CBETanTRGSweep(rho,H,Para,env);
end
info.Mode=mode; info.SweepSeconds=toc(timer);
info.MaxBond=max(cellfun(@(a)size(a,1),rho.A));
end
