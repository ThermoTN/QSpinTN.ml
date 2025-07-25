color_list = {[0.9290 0.6940 0.1250], [0.8500 0.3250 0.0980], [0.6350 0.0780 0.1840], [0.3010 0.7450 0.9330], [0 0.4470 0.7410]};
Rslt = loadRslt(50, 'tanTRG');
[Fe, En, Cm, S] = FreeFermion_1Ddst(50, Rslt.beta());
loglog(Rslt.beta, abs(Rslt.En - En'), '-o', 'linewidth', 2, 'color', color_list{1}); hold on

Rslt = loadRslt(100, 'tanTRG');
[Fe, En, Cm, S] = FreeFermion_1Ddst(50, Rslt.beta());
loglog(Rslt.beta, abs(Rslt.En - En'), '-o', 'linewidth', 2, 'color', color_list{2}); hold on

Rslt = loadRslt(200, 'tanTRG');
[Fe, En, Cm, S] = FreeFermion_1Ddst(50, Rslt.beta());
loglog(Rslt.beta, abs(Rslt.En - En'), '-o', 'linewidth', 2, 'color', color_list{3}); hold on

Rslt = loadRslt(50, 'XTRG');
[Fe, En, Cm, S] = FreeFermion_1Ddst(50, Rslt.beta());
loglog(Rslt.beta, abs(Rslt.En - En), '-s', 'linewidth', 2, 'color', color_list{4}); hold on

Rslt = loadRslt(100, 'XTRG');
[Fe, En, Cm, S] = FreeFermion_1Ddst(50, Rslt.beta());
loglog(Rslt.beta, abs(Rslt.En - En), '-s', 'linewidth', 2, 'color', color_list{5}); hold on


xlabel('$\beta$', 'Interpreter', 'latex')
ylabel('$\delta e$', 'Interpreter', 'latex')
set(gca, 'linewidth', 1.5, 'fontsize', 20, 'XColor', 'k', 'YColor', 'k', 'fontname', 'times new roman')
legend({'$D=50$', '$D=100$'}, 'Interpreter', 'latex', 'location', 'northwest', 'box', 'off')
legend({'tanTRG $D=50$', '$D=100$', '$D=200$', 'XTRG $D=50$', '$D=100$'}, 'Interpreter', 'latex', 'location', 'northwest', 'box', 'off')
ylim([1e-16, 1e-3])
hold off

function Rslt = loadRslt(D, mbsolver)
switch mbsolver
    case 'tanTRG'
        load(['Rslt/Model=XXZtest-Jxy=1-Jz=0-h=[0 0 0]-L=50-BC=OBC-D=', num2str(D), '.mat'])
        Rslt.beta = Rslt.beta(2:end);
        Rslt.En = Rslt.En(2:end);
        Rslt.LnZ = Rslt.LnZ(2:end);
        Rslt.EE = Rslt.EE(2:end);
        Rslt.M = Rslt.M(2:end);
        Rslt.TrunErr = Rslt.TrunErr(2:end);
    case 'XTRG'
        load(['Rslt/Model=XXZtest-Jxy=1-Jz=0-h=[0 0 0]-L=50-BC=OBC-D=', num2str(D), '-XTRG.mat'])
end
end
function [Fe, En, Cm, S] = FreeFermion_1Ddst(L, beta)
% exact result of OBC*OBC free spinless fermions
% 
% if strcmp(varargin, 'PBC')
%    epsilon = -2*cos((1:L)*pi/(L+1)) - 2*cos((0:W-1).'*2*pi/W);
% else
%     epsilon = -2*cos((1:L)*pi/(L+1)) - 2*cos((1:W).'*pi/(W+1));
% end

epsilon = -cos((1:L)*pi/(L+1));

epsilon = reshape(epsilon, [numel(epsilon), 1]);

[Fe, En, Cm, S] = calThermo(epsilon, beta);

end

function [Fe, En, Cm, S] = calThermo(epsilon, beta) 

% calculate thermodynamics quantities from spectrum of fermions


epsilon = reshape(epsilon, [numel(epsilon), 1]); % column vector
L = numel(epsilon);

beta = reshape(beta, [1, numel(beta)]); % row vector

% Fe = -ln2/beta - 1/beta*L sum(lncosh(beta*epsilon_k/2))
Fe = - (sum(my_lncosh(epsilon*beta/2), 1)/L + log(2))./beta;
Fe = Fe.';

% En = - 1/2L sum(epsilon_k tanh(beta*epsilon_k/2))
En = - sum(diag(epsilon)*tanh(epsilon*beta/2), 1)/(2*L);
En = En.';

% Cm = beta^2/2L sum(epsilon_k^2/(1 + cosh(beta*epsilon_k) ))
Cm = sum(epsilon.^2./(1 + cosh(epsilon*beta)), 1).*beta.^2/(2*L);
Cm = Cm';

% S = beta*(En - Fe)
S = beta'.*(En - Fe);

end

function [y] =  my_lncosh(x)

% fixed bugs that cosh(x) is too large to output NaN
y = log(cosh(x));

a = isinf(y);

% x is large enough, lncosh(x) ~ x - ln2
y(a) = x(a) - log(2);

end