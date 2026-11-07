function sweepTol = GetDMRGSweepTol(Para)
% Other models and archived configurations may predate DMRGSweepTol.
sweepTol = 1e-8;
if isfield(Para,'DMRGSweepTol')
    sweepTol = Para.DMRGSweepTol;
end
validateattributes(sweepTol,{'double','single'},{'scalar','real','finite','positive'}, ...
    mfilename,'Para.DMRGSweepTol');
end
