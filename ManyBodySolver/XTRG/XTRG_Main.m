function [ Rslt ] = XTRG_Main( Para, H, Id, Op, ST )
taup = Para.tau0 * 2 ^ (ST/4);
Para.tau0 = taup;
fileID = Para.fileID;
[rho_init, Op] = SETTN(Para, H, Id, Op, fileID);
[~, Rslt] = ExpEvo_conj(Para, rho_init, Op, fileID);
end

