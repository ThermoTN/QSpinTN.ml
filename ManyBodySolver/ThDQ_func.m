function [Cm, beta] = ThDQ_func(Rslt)
skip = 4;

beta_list = 0.0025.*2.^(0:1:15).*2^0.1;
for int = 0.3:0.2:0.9
    beta_list = [beta_list, 0.0025.*2.^(0:1:15).*2^int];
end
beta_list = sort(beta_list);
loc1 = find(beta_list > Rslt.beta(skip+1), 1);
loc2 = find(beta_list > Rslt.beta(end-skip), 1);
if isempty(loc2)
    loc2 = length(beta_list);
end
beta = beta_list(loc1:loc2-1);
En = real(Rslt.En);
pp = csaps(log(Rslt.beta), En); % csaps is better than spline!
Cm = zeros([length(beta), 1]);
for it = 1:length(beta)
    Cm(it) = getCm(pp, beta(it));
end
end

function Cm = getCm(pp, beta)
lgbeta = log(beta);

loc = find(log(beta) >= pp.breaks);
loc = loc(end);

Cm = -beta * (pp.coefs(loc, 1) * 3 * (lgbeta - pp.breaks(loc)) ^ 2 ...
    + pp.coefs(loc, 2) * 2 * (lgbeta - pp.breaks(loc)) ...
    + pp.coefs(loc, 3));
end