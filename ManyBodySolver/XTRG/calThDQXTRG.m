function [ ThDQ ] = calThDQXTRG( Para, rho, Op )
% function [ ThDQ ] = calThDQ( Para, rho, H )

if norm(Para.Field.h) ~= 0
    ThDQ.M = real(InnerProd(rho, Op.Sm));
else
    ThDQ.M = 0;
end

ThDQ.En = real(InnerProd(rho, Op.H));

end

