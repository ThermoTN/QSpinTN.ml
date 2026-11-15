function opts = GetCBEOptions(Para)
% Pre-expansion target follows tanTRG_Spin/BondExpansionfunc/PreBE_RSVD.
opts=struct('Ratio',.1,'Oversampling',8,'PowerIterations',1,'SvdTol',1e-12);
names=fieldnames(opts);
for i=1:numel(names)
    field=['CBE',names{i}];
    if isfield(Para,field), opts.(names{i})=Para.(field); end
end
validateattributes(opts.Ratio,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(opts.Oversampling,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(opts.PowerIterations,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(opts.SvdTol,{'numeric'},{'scalar','nonnegative','<',1,'finite'});
opts.Cap=ceil((1+opts.Ratio)*Para.MCrit);
if isfield(Para,'CBEExpansionCap') && ~isempty(Para.CBEExpansionCap)
    validateattributes(Para.CBEExpansionCap,{'numeric'},{'scalar','integer','>=',Para.MCrit,'finite'});
    opts.Cap=min(opts.Cap,Para.CBEExpansionCap);
end
end
