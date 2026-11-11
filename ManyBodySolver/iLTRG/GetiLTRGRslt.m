function Rslt = GetiLTRGRslt(Para, ThDQ)
%GETILTRGRSLT Natural-unit result on one consistent, non-extrapolated grid.
% M is bare <S> along the effective field; use an explicit unit conversion.
if ~ismember(ThDQ,{'Cm','Chi','Cm&Chi'})
    error('qspintn:iltrg:Request','ThDQ must be Cm, Chi or Cm&Chi.');
end
[~,~,raw]=iTEBD(Para);
if strcmp(ThDQ,'Chi')
    indices=(1:Para.N_max)';
else
    if Para.N_max<3, error('qspintn:iltrg:Request','Cm needs at least three steps.'); end
    indices=(2:Para.N_max-1)';
end
Rslt=struct('T',raw.T_K(indices),'beta',raw.beta(indices), ...
    'LnZ',raw.lnZ(indices),'M',raw.spin_projection(indices), ...
    'Cm',raw.C_natural(indices));
Rslt.grid='native; no interpolation or extrapolation';
Rslt.units='natural per spin; T is in the model energy unit';
end
