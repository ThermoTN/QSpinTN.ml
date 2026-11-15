classdef DMRGEnvironmentStore < SweepEnvironmentStore
    % Preserve the public DMRG store API and parameter names.
    methods
        function obj = DMRGEnvironmentStore(Para)
            obj@SweepEnvironmentStore(Para,'DMRG');
        end
    end
end
