maxNumCompThreads(4);
load('YMGO_best7_DMRG.mat');
[Spin] = Observe_ezLv(RsltChi{1,1}{1,1}.T.A);
save Spin_0T Spin
load('YMGO_DMRG_1T.mat');
[Spin] = Observe_ezLv(RsltChi{1,1}{1,1}.T.A);
save Spin_1T Spin
load('YMGO_DMRG_2T.mat');
[Spin] = Observe_ezLv(RsltChi{1,1}{1,1}.T.A);
save Spin_2T Spin
load('YMGO_DMRG_3T.mat');
[Spin] = Observe_ezLv(RsltChi{1,1}{1,1}.T.A);
save Spin_3T Spin
load('YMGO_DMRG_4T.mat');
[Spin] = Observe_ezLv(RsltChi{1,1}{1,1}.T.A);
save Spin_4T Spin