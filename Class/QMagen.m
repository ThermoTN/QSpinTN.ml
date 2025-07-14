classdef QMagen < matlab.mixin.CustomDisplay
    % A data type for QMagen Optimaization with properties
    % [.] has default values or can be computed from other variables.
    %
    % .Config
    %       .ManyBodySolver
    %       .ModelName
    %       .Mode
    %
    % .ModelConf
    %       .[ModelName]
    %       .[ModelName_Full]
    %       .[IntrcMap]
    %       .[LocalSpin]
    %       .[Para_Name]
    %       .[Para_EnScale]
    %       .Para_Range         % For 'OPT' & 'LOSS'
    %       .[gFactor_Num]
    %       .[gFactor_Type]
    %       .[gFactor_Name]
    %       .[gFactor_Vec]
    %       .gFactor_Range      % For 'OPT' & 'LOSS'
    %
    % .Lattice
    %       .[L]
    %       .Lx
    %       .Ly
    %       .BCX
    %       .BCY
    %
    % .Field
    %       .B
    %       .h
    %
    % .LossConf
    %       .WeightList
    %       .Type
    %       .Design
    %
    % .Setting
    %       .PLOTFLAG
    %       .EVOFLAG
    %       .SAVEFLAG
    %       .SAVENAME
    % -------------------------------------------------
    % Usage:
    %       QMagenConf(Conf, ModelConf, Lattice, LossConf, Setting, Restart, 'Thermo data name', ThermoData)
    %       QMagenConf(Conf, ModelConf, Lattice)
    % -------------------------------------------------
    % BBC, YG@BUAA, Mar16
    
    properties % (Access = private)
        Config
        ModelConf
        ModelParaValue
        Lattice
        Field
        LossConf
        Setting
        Restart = {};
        CmData = {};
        ChiData = {};
    end
    
    methods
        function obj = QMagen(Config, ModelConf, Lattice, varargin)
            if nargin == 9 || nargin == 11
                obj.Config = Config;
                obj.ModelConf = ModelConf;
                obj.Lattice = Lattice;
                obj.LossConf = varargin{1};
                obj.Setting = varargin{2};
                len = length(varargin);
                % //whether the many-body solver is available for the model or not
                if ~ismember(Config.ManyBodySolver, obj.ModelConf.AvlbSolver)
%                    warning('Para.ManyBodySolver not within ModelConf.AvlbSolver!');
%                    pause;
                end
                if len ~= 4 && len ~= 6 && len ~= 8
                    error('Illegal thermodynamic data import format!')
                else
                    for i = 1:1:(len-2)/2
                        switch varargin{2 * i + 1}
                            case {'C', 'c', 'Cm', 'cm'}
                                obj.CmData = varargin{2 * i + 2};
                            case {'Chi', 'chi'}
                                obj.ChiData = varargin{2 * i + 2};
                            case {'Restart'}
                                obj.Restart = varargin{2 * i + 2};
                            otherwise
                                error('Illegal data import format!')
                        end
                    end
                end
            elseif nargin == 4
                obj.Config = Config;
                obj.ModelConf = ModelConf;
                obj.Lattice = Lattice;
                obj.Field = varargin{1};
            elseif nargin == 8 || nargin == 6
                obj.Config = Config;
                obj.ModelConf = ModelConf;
                obj.Lattice = Lattice;
                obj.Setting = varargin{1};
                len = length(varargin);
                for i = 1:1:(len-1)/2
                    switch varargin{2 * i}
                        case {'C', 'c', 'Cm', 'cm'}
                            obj.CmData = varargin{2 * i + 1};
                        case {'Chi', 'chi'}
                            obj.ChiData = varargin{2 * i + 1};
                        otherwise
                            error('Illegal data import format!')
                    end
                end
            else
                error('Illegal import format!')
            end
        end
        
    end
    
    methods (Access=protected)
        function displayScalarObject(obj)
            if isempty(obj.Config)
                fprintf('(empty ThermoData)\n');
            else
                switch obj.Config.Mode
                    case 'OPT'
                        fprintf('\nQMagen Configuration:\n');
                    case 'ThDQ'
                        fprintf('\nQSpinLib Configuration:\n');
                end
                fprintf('-------------------------------------------------------------------------------------------\n')
                fprintf('\tModel Name           |  %s\n', obj.ModelConf.ModelName_Full);
                fprintf('\tMany-Body Solver     |  %s\n', obj.Config.ManyBodySolver);
                if strcmp(obj.Config.Mode, 'OPT')
                    fprintf('\tOptimizer            |  %s\n', 'Bayesian Optimizer');
                end
                if strcmp(obj.Config.Mode, 'ThDQ')
                    fprintf('\tWorking mode         |  %s\n', 'Finite-temperature simulation');
                end
                fprintf('\tStore folder         |  %s\n', ['../Tmp/tmp_', obj.Config.TStr, '/']);
                fprintf('-------------------------------------------------------------------------------------------\n')
                fprintf('Lattice Geometry:\n');
                fprintf('\tL                    |  %d\n', obj.Lattice.L);
                if isfield(obj.Lattice, 'Lx')
                    fprintf('\tLx                   |  %d\n', obj.Lattice.Lx);
                end
                if isfield(obj.Lattice, 'Ly')
                    fprintf('\tLy                   |  %d\n', obj.Lattice.Ly);
                end
                if isfield(obj.Lattice, 'Lx')
                    fprintf('\tBCX                  |  %s\n', obj.Lattice.BCX);
                end
                if isfield(obj.Lattice, 'Lx')
                    fprintf('\tBCY                  |  %s\n', obj.Lattice.BCY);
                end
                if isfield(obj.Lattice, 'BC')
                    fprintf('\tBC                   |  %s\n', obj.Lattice.BC);
                end
                fprintf('-------------------------------------------------------------------------------------------\n')
                switch obj.Config.Mode
                    case 'OPT'
                        fprintf('Parameters Optimization Range:\n');
                        for i = 1:1:length(obj.ModelConf.Para_Name)
                            fprintf('\t%-18s   |  ', obj.ModelConf.Para_Name{i})
                            if ischar(obj.ModelConf.Para_Range{i})
                                fprintf('%s\n', obj.ModelConf.Para_Range{i})
                            elseif length(obj.ModelConf.Para_Range{i}) == 2
                                fprintf('[%.3f, %.3f]\n', obj.ModelConf.Para_Range{i})
                            elseif length(obj.ModelConf.Para_Range{i}) == 1
                                fprintf('%.3f\n', obj.ModelConf.Para_Range{i})
                            else
                                error('??????????\n')
                            end
                        end
                        fprintf('-------------------------------------------------------------------------------------------\n')
                    case {'LOSS', 'CALC-Cm', 'CALC-Chi'}
                        fprintf('Parameters Value:\n');
                        for i = 1:1:length(obj.ModelConf.Para_Name)
                            fprintf('\t%-18s   |  %.3f\n', obj.ModelConf.Para_Name{i}, obj.ModelConf.Para_Value{i})
                        end
                        fprintf('-------------------------------------------------------------------------------------------\n')
                        % fprintf('Parameters Value:\n');
                        fprintf('\t%-18s   |  [%.3f, %.3f, %.3f] \n', 'Field (Tesla)', obj.Field.B)
                        fprintf('\t%s |  [%.3f, %.3f, %.3f] \n', 'Field (Natural Unit)', obj.Field.h)
                        try
                            fprintf('\t%-18s   |  [%.3f, .%.3f, %.3f] \n', 'g', obj.ModelParaValue.g)
                        catch
                        end
                        fprintf('-------------------------------------------------------------------------------------------\n')
                    case {'ThDQ'}
                        fprintf('Parameters Value: \n');
                        for i = 1:1:length(obj.ModelConf.Para_Name)
                            if strcmp(obj.ModelConf.Para_Unit{i}, 'K')
                                unit = 'Kelvin';
                                CS = eval(['obj.ModelParaValue.', obj.ModelConf.Para_Name{i}]) * obj.ModelParaValue.ES;
                            else
                                unit = ['[', obj.ModelConf.Para_EnScale, ']'];
                                CS = eval(['obj.ModelParaValue.', obj.ModelConf.Para_Name{i}]);
                            end
                            fprintf('\t%-18s   |  %.3f %s\n', obj.ModelConf.Para_Name{i}, CS, unit)
                        end
                        fprintf('\t%-18s   |  [%.3f, %.3f, %.3f] Tesla\n', 'Field', obj.Field.B)
                        fprintf('-------------------------------------------------------------------------------------------\n')
                        fprintf('Parameters Value: (Naturel Unit)\n');
                        for i = 1:1:length(obj.ModelConf.Para_Name)
                            fprintf('\t%-18s   |  %.3f\n', obj.ModelConf.Para_Name{i}, eval(['obj.ModelParaValue.', obj.ModelConf.Para_Name{i}]))
                        end
                        % fprintf('Parameters Value:\n');
                        fprintf('\t%-18s   |  %.3f \n', 'Field', obj.Field.h)
%                         try
%                             fprintf('\t%-18s   |  [%.3f, %.3f, %.3f] \n', 'g', obj.ModelParaValue.g)
%                         catch
%                         end
                        fprintf('-------------------------------------------------------------------------------------------\n')
                end
            end
        end
        
        function displayNonScalarObject(objAry)
            for obj=objAry
                displayScalarObject(obj);
            end
        end
        
    end
    
    
end
