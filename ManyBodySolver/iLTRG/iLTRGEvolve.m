function [raw, state, info] = iLTRGEvolve(Para, Hab, Hba, options)
%ILTRGEVOLVE Shared two-site iLTRG evolution with optional durable state.
% Para.tau is dimensionless; Para.energy_scale_K defaults to 1. Output beta
% and T_K use the corresponding physical scale. No project files are made.
% Copyright (c) 2025 ThermoTN. MIT License; see repository LICENSE.
if nargin<4, options=struct(); end
options = validate_options(options);
if ~isfield(Para,'energy_scale_K'), Para.energy_scale_K=1; end
if ~isfield(Para,'transfer_tolerance'), Para.transfer_tolerance=1e-12; end
if ~isfield(Para,'transfer_max_iterations'), Para.transfer_max_iterations=10000; end
for name={'tau','energy_scale_K','transfer_tolerance'}
    validateattributes(Para.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
for name={'d','D_max','N_max','transfer_max_iterations'}
    validateattributes(Para.(name{1}),{'numeric'},{'scalar','real','finite','positive','integer'});
end
validateattributes(Para.Field.h,{'numeric'},{'vector','numel',3,'real','finite'});
for h={Hab,Hba}
    matrix=h{1};
    if ~isequal(size(matrix),[Para.d^2,Para.d^2]) || any(~isfinite(matrix(:))) || ...
            norm(matrix-matrix','fro')>1e-12*max(1,norm(matrix,'fro'))
        error('qspintn:iltrg:Hamiltonian','Expected a finite Hermitian d^2 by d^2 Hamiltonian.');
    end
end
sources = iLTRGSourceIdentity();
science = struct('schema_version','qspintn.iltrg.state.v1','parameters',Para, ...
    'Hab_real',real(Hab),'Hab_imag',imag(Hab),'Hba_real',real(Hba),'Hba_imag',imag(Hba), ...
    'sources',sources);
identity=hash_bytes(uint8(unicode2native(jsonencode(science),'UTF-8')));
state=options.state;
if ~isempty(state) && ~isempty(options.checkpoint_path) && isfile(options.checkpoint_path)
    error('qspintn:iltrg:StateConflict','Choose an in-memory state or a saved checkpoint, not both.');
end
if isempty(state) && ~isempty(options.checkpoint_path) && isfile(options.checkpoint_path)
    saved=load(options.checkpoint_path,'checkpoint');
    if ~isfield(saved,'checkpoint'), error('qspintn:iltrg:CheckpointMismatch','Missing checkpoint state.'); end
    state=saved.checkpoint;
end
iteration=0; lgnorm=0; elapsed_prior=0; env=struct();
raw=empty_raw(Para.N_max);
if ~isempty(state)
    required={'identity','Ta','Tb','La','Lb','iteration','lgnorm','raw','env','elapsed_seconds'};
    if ~isstruct(state) || ~isscalar(state) || ~all(isfield(state,required)) || ...
            ~strcmp(state.identity,identity)
        error('qspintn:iltrg:CheckpointMismatch','State model, numerical settings, units or numerical sources differ.');
    end
    validateattributes(state.iteration,{'numeric'},{'scalar','integer','nonnegative','<=',Para.N_max});
    Ta=state.Ta; Tb=state.Tb; La=state.La; Lb=state.Lb;
    iteration=state.iteration; lgnorm=state.lgnorm;
    raw=state.raw; env=state.env; elapsed_prior=state.elapsed_seconds;
else
    [Ta,Tb,La,Lb]=InitId(Para);
end
resumed_from_iteration=iteration;
[Uab,Uba]=EvoGate(Hab,Hba,Para);
physical_tau=Para.tau/Para.energy_scale_K;
need_m=norm(Para.Field.h)>0;
status='running'; t0=tic;
notify_progress();
try
    for it=iteration+1:Para.N_max
        if ~isempty(options.should_stop) && options.should_stop()
            status='paused'; break;
        end
        fprintf(options.log_fid,'It: %6d| beta_K_inv: %.12g| ',it,2*physical_tau*it);
        [Ta,Lb,Tb,Ns]=EvoRG(La,Ta,Lb,Tb,Uab,Para,options.log_fid);
        lgnorm=lgnorm+log(positive_real(Ns,'first normalization'));
        [Tb,La,Ta,Ns]=EvoRG(Lb,Tb,La,Ta,Uba,Para,options.log_fid);
        lgnorm=lgnorm+log(positive_real(Ns,'second normalization'));
        if any(~isfinite(Ta(:))) || any(~isfinite(Tb(:))) || ...
                any(~isfinite(La(:))) || any(~isfinite(Lb(:)))
            error('qspintn:iltrg:NonFinite','Non-finite evolved tensor at step %d.',it);
        end
        [lambda,m,diagnostic,env]=iLTRGObserve(Ta,La,Tb,Lb,Para,need_m,env);
        raw.beta(it)=2*physical_tau*it;
        raw.T_K(it)=1/raw.beta(it);
        raw.lnZ(it)=checked_real(lgnorm+log(positive_real(lambda,'transfer eigenvalue'))/2,'lnZ');
        raw.spin_projection(it)=checked_real(m,'spin projection');
        raw.actual_D_a(it)=size(La,1); raw.actual_D_b(it)=size(Lb,1);
        raw.transfer{it}=diagnostic;
        if it>=3
            j=it-1;
            raw.C_natural(j)=raw.beta(j)^2* ...
                (raw.lnZ(j-1)+raw.lnZ(j+1)-2*raw.lnZ(j))/(4*physical_tau^2);
        end
        iteration=it;
        fprintf(options.log_fid,'LnZ: %.16g| M_spin: %.16g| Da: %d| Db: %d\n', ...
            raw.lnZ(it),m,raw.actual_D_a(it),raw.actual_D_b(it));
        notify_progress();
        if mod(iteration,options.checkpoint_interval)==0, save_checkpoint(); end
        if iteration-resumed_from_iteration>=options.pause_after_steps && iteration<Para.N_max
            status='paused'; break;
        end
    end
    if iteration==Para.N_max, status='completed'; end
    if strcmp(status,'running'), status='paused'; end
    save_checkpoint();
    notify_progress();
catch err
    % Evolving a bond is not an atomic step. Preserve the last durable safe
    % checkpoint rather than saving partially advanced tensors on failure.
    status='failed';
    try, notify_progress(); catch, end
    rethrow(err);
end
info=progress(); info.sources=sources; info.science=science;
state=make_state();
names=fieldnames(raw);
for k=1:numel(names), raw.(names{k})=raw.(names{k})(1:iteration); end

    function value=make_state()
        value=struct('schema_version','qspintn.iltrg.state.v1','identity',identity, ...
            'science',science,'Ta',Ta,'Tb',Tb,'La',La,'Lb',Lb, ...
            'iteration',iteration,'lgnorm',lgnorm,'raw',raw,'env',env, ...
            'elapsed_seconds',elapsed_prior+toc(t0));
    end
    function save_checkpoint()
        if isempty(options.checkpoint_path), return; end
        checkpoint=make_state(); %#ok<NASGU>
        parent=fileparts(options.checkpoint_path);
        tmp=[tempname(parent),'.mat'];
        cleanup=onCleanup(@() remove_file(tmp)); %#ok<NASGU>
        save(tmp,'checkpoint','-v7.3');
        [ok,msg]=movefile(tmp,options.checkpoint_path,'f');
        if ~ok, error('qspintn:iltrg:CheckpointWrite','%s',msg); end
    end
    function p=progress()
        p=struct('status',status,'iteration',iteration,'N_max',Para.N_max, ...
            'identity',identity,'resumed_from_iteration',resumed_from_iteration, ...
            'D_max',Para.D_max,'elapsed_seconds',elapsed_prior+toc(t0));
        if iteration>0
            p.T_K=raw.T_K(iteration); p.actual_D_a=raw.actual_D_a(iteration);
            p.actual_D_b=raw.actual_D_b(iteration); p.transfer=raw.transfer{iteration};
        end
    end
    function notify_progress()
        if ~isempty(options.on_progress), options.on_progress(progress()); end
    end
end

function opts=validate_options(opts)
defaults=struct('checkpoint_path','','checkpoint_interval',10,'pause_after_steps',inf, ...
    'log_fid',1,'on_progress',[],'should_stop',[],'state',[]);
if ~isstruct(opts) || ~isscalar(opts) || ~all(ismember(fieldnames(opts),fieldnames(defaults)))
    error('qspintn:iltrg:Options','Unsupported numerical runtime options.');
end
names=fieldnames(defaults);
for k=1:numel(names)
    if ~isfield(opts,names{k}), opts.(names{k})=defaults.(names{k}); end
end
validateattributes(opts.checkpoint_interval,{'numeric'},{'scalar','real','finite','integer','positive'});
if opts.pause_after_steps~=inf
    validateattributes(opts.pause_after_steps,{'numeric'},{'scalar','real','finite','integer','positive'});
end
validateattributes(opts.log_fid,{'numeric'},{'scalar','integer','positive'});
for name={'on_progress','should_stop'}
    if ~isempty(opts.(name{1})) && ~isa(opts.(name{1}),'function_handle')
        error('qspintn:iltrg:Options','%s must be a function handle.',name{1});
    end
end
if isstring(opts.checkpoint_path) && isscalar(opts.checkpoint_path), opts.checkpoint_path=char(opts.checkpoint_path); end
if ~ischar(opts.checkpoint_path), error('qspintn:iltrg:Options','checkpoint_path must be text.'); end
if ~isempty(opts.checkpoint_path)
    if ~java.io.File(opts.checkpoint_path).isAbsolute()
        error('qspintn:iltrg:Options','checkpoint_path must be absolute.');
    end
    if ~isfolder(fileparts(opts.checkpoint_path))
        error('qspintn:iltrg:Options','The caller must create the checkpoint directory.');
    end
end
end

function raw=empty_raw(n)
raw=struct();
for name={'beta','T_K','lnZ','spin_projection','C_natural','actual_D_a','actual_D_b'}
    raw.(name{1})=nan(n,1);
end
raw.transfer=cell(n,1);
end
function value=checked_real(value,label)
if ~isscalar(value) || ~isfinite(value) || abs(imag(value))>1e-11*max(1,abs(real(value)))
    error('qspintn:iltrg:Numerics','Nonfinite or materially complex %s.',label);
end
value=real(value);
end
function value=positive_real(value,label)
value=checked_real(value,label);
if value<=0, error('qspintn:iltrg:Normalization','Nonpositive %s.',label); end
end
function digest=hash_bytes(bytes)
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(typecast(uint8(bytes(:)),'int8'));
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function remove_file(location)
if isfile(location), delete(location); end
end
