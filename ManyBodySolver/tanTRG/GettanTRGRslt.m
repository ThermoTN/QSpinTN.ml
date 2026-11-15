function Rslt = GettanTRGRslt(Para)
% SETTN initialization followed by selectable thermal sweeps.
[H,Id,Op] = InitHam(Para);
[rho,Op] = SETTN(Para,H,Id,Op,Para.fileID);
Rslt = tanTRGMain(Para,H,rho,Op);
end
