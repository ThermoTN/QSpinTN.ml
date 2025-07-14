function [ Rslt ] = caltdvpThDQ( Para, rho, Op, Rslt )

if norm(Para.Field.h) ~= 0
    Rslt.M(end + 1) = real(InnerProd(rho, Op.Sm))/Para.L;
else
    Rslt.M(end + 1) = 0;
end

end

