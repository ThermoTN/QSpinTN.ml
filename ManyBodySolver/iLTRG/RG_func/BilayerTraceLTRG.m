function maxeig = BilayerTraceLTRG(Ta, La, Tb, Lb)
%BILAYERTRACELTRG Two-site transfer eigenvalue from the shared fixed point.
[maxeig, ~] = iLTRGObserve(Ta,La,Tb,Lb,struct(),false);
end
