function energy = GetEn(EnV,H,rho)
% Bare normalized energy, without exp(H.lgnorm) or per-site division.
d=size(rho.A{1},2); a=reshape(rho.A{1},1,[],d*d);
w=ThermalSweepLayout.lift(reshape(H.A{1},1,[],d,d),d);
apply=TensorSweep.operator(1,w,permute(EnV,[1,3,2]),size(a));
y=apply(a); energy=a(:)'*y(:);
end
