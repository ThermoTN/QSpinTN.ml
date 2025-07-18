function [ En ] = Observe( MPS, Ham, Para )

% //DMRG runtime parameters
Para.Ver = 'Memory';
% max iterations of MPO varitional product
Para.VariProd_step_max = 10000;
% max iterations of MPO varitional sum
Para.VariSum_step_max = 10000;
% max expensian order of SETTN
Para.SETTN_init_step_max = 10000;


[ C, ns ] = VariProdMPS( Para, MPS.A, Ham.A );
res = real(InnerProdMPS(C, MPS.A));
% if ~isempty(varargin)
%     lognorm = lognorm + 2 * varargin{1};
% end
En = ns * res * exp(Ham.lgnorm);

end

