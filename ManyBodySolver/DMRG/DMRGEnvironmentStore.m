classdef DMRGEnvironmentStore < handle
    % Rolling sweep environments. A slot holds either its left or right block.
    properties (SetAccess=private)
        Directory = ''
    end
    properties (Access=private)
        Blocks
        Sides
        Keep = true
    end
    methods
        function obj = DMRGEnvironmentStore(Para)
            mode = 'memory';
            if isfield(Para,'DMRGEnvStorage'), mode = Para.DMRGEnvStorage; end
            mode = validatestring(mode,{'memory','disk'});
            obj.Sides = repmat(' ',1,Para.L);
            if strcmp(mode,'memory')
                obj.Blocks = cell(1,Para.L); return
            end
            if isfield(Para,'DMRGKeepEnv')
                validateattributes(Para.DMRGKeepEnv,{'logical'},{'scalar'});
                obj.Keep = Para.DMRGKeepEnv;
            end
            root = fullfile(tempdir,'qspintn_dmrg');
            if isfield(Para,'DMRGEnvRoot')
                root = Para.DMRGEnvRoot;
                if isstring(root) && isscalar(root), root = char(root); end
                validateattributes(root,{'char'},{'row','nonempty'},mfilename,'Para.DMRGEnvRoot');
            end
            [ok,msg] = mkdir(root); assert(ok,'DMRG:EnvDirectory','%s: %s',root,msg);
            [ok,attr] = fileattrib(root);
            assert(ok,'DMRG:EnvDirectory','Cannot resolve environment root %s.',root);
            root = attr.Name;
            % mkdir on a fresh path also excludes accidental reuse by concurrent jobs.
            [~,token] = fileparts(tempname(root));
            job = regexprep(getenv('SLURM_JOB_ID'),'[^a-zA-Z0-9_-]','_');
            if isempty(job), job = 'local'; end
            path = fullfile(root,['dmrg_',job,'_',token]);
            [ok,msg,id] = mkdir(path);
            assert(ok && isempty(id),'DMRG:EnvDirectory','Cannot create fresh directory %s: %s',path,msg);
            obj.Directory = path;
            fprintf('DMRG environment directory: %s\n',path);
        end
        function put(obj,slot,side,block)
            if isempty(obj.Directory)
                obj.Blocks{slot} = block;
            else
                % Commit a whole block before replacing the previous direction.
                file = obj.filename(slot); tmp = [file,'.partial.mat'];
                save(tmp,'block','slot','side','-v7.3','-nocompression');
                [ok,msg] = movefile(tmp,file,'f');
                assert(ok,'DMRG:EnvWrite','Cannot replace %s: %s',file,msg);
            end
            obj.Sides(slot) = side;
        end
        function block = get(obj,slot,side)
            assert(obj.Sides(slot)==side,'DMRG:EnvSide','Stale or missing %s environment at %d.',side,slot);
            if isempty(obj.Directory)
                block = obj.Blocks{slot};
            else
                entry = load(obj.filename(slot),'block','side','slot');
                assert(entry.side==side && entry.slot==slot,'DMRG:EnvFile','Unexpected environment at %d.',slot);
                block = entry.block;
            end
        end
        function reset(obj)
            % Fused two-site environments and bare single-site environments
            % have different shapes, even at the same slot and sweep side.
            if isempty(obj.Directory)
                obj.Blocks = cell(size(obj.Sides));
            else
                for slot=find(obj.Sides~=' ')
                    delete(obj.filename(slot));
                end
            end
            obj.Sides(:) = ' ';
        end
        function finish(obj)
            if ~obj.Keep && ~isempty(obj.Directory) && isfolder(obj.Directory)
                [ok,msg] = rmdir(obj.Directory,'s');
                if ~ok, warning('DMRG:EnvCleanup','%s: %s',obj.Directory,msg); end
            end
        end
    end
    methods (Access=private)
        function path = filename(obj,slot)
            path = fullfile(obj.Directory,sprintf('env_%04d.mat',slot));
        end
    end
end
