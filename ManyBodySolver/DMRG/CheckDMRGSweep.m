function [converged,deltaS,deltaE] = CheckDMRGSweep(energy,previous,entropy,previousEntropy,tol)
% Shared historical two-site criterion: relative E AND change of max bond S.
deltaE = real(energy)-real(previous);
deltaS = abs(max(entropy(:))-max(previousEntropy(:)));
converged = isfinite(previous) && ...
    abs(deltaE)<tol*max(eps,abs(previous)) && deltaS<tol;
end
