% load 'HeiModel_h=[0 0 0].mat'
% 
% 
% data(:,1) = 1./Rslt.beta(1:end-1);
% data(:,2) = Rslt.Cm(1:end-1) * 8.314;
% 
% plot(data(:,1), data(:,2))
% 
% save('XXZChain_Jxy=1_Jz=1.5_g=2_Cm_0T.mat', 'data')


% load 'HeiModel_h=[0 0 0.134362625805593].mat'
% 
% data(:,1) = 1./Rslt.beta(1:end-1);
% data(:,2) = Rslt.M(1:end-1) * 1.178707921721845 * 4;
% 
% save('XXZChain_Jxy=1_Jz=1.5_g=2_Chi_Bz=0.1T.mat', 'data')
% 
% load 'HeiModel_h=[0.134362625805593 0 0].mat'
% 
% data(:,1) = 1./Rslt.beta(1:end-1);
% data(:,2) = Rslt.M(1:end-1) * 1.178707921721845 * 4;
% 
% save('XXZChain_Jxy=1_Jz=1.5_g=2_Chi_Bx=0.1T.mat', 'data')


load 'HeiModel_h=[0 0 0.1].mat'
semilogx(1./Rslt.beta, Rslt.M/0.1)

xlim([0, 10])
xlabel('T')
ylabel('M/h')