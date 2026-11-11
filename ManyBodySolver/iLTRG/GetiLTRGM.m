function M = GetiLTRGM(Ta, La, Tb, Lb, Para)
%GETILTRGM Natural spin projection from the shared converged observable.
[~, M] = iLTRGObserve(Ta,La,Tb,Lb,Para,true);
end
