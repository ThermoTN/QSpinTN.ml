function [ Para ] = PassPara( QMagenConf, K_min )
% function [Para] = GetPara(Model, Field, K_min)
% Set the parameter of the problem.

% Para.tau0 = QMagenConf.Config.tau0;
Para.ManyBodySolver = QMagenConf.Config.ManyBodySolver; 
Para.TStr = QMagenConf.Config.TStr;
Para.TStr_log = QMagenConf.Config.TStr_log;

% //pass interaction map (ED, XTRG) or Trotter gates (iLTRG) to solvers
try
    Para.IntrcMap_Name = QMagenConf.ModelConf.IntrcMap;
end
try
    Para.Trotter_Name = QMagenConf.ModelConf.Trotter;
end

Para.d = eval(QMagenConf.ModelConf.LocalSpin) * 2 + 1;
Para.L = QMagenConf.Lattice.L;
Para.Geo = QMagenConf.Lattice;

% if strcmp(QMagenConf.Config.Mode, 'ThDQ') && strcmp(QMagenConf.Config.ManyBodySolver, 'XTRG')
%     Para.Parallel = true;
% else
%     Para.Parallel = false;
% end
% ====================================================
% switch QMagenConf.Config.Mode
%     case {'OPT', 'LOSS', 'ThDQ'}
%         global PLOTFLAG
%         try
%             PLOTFLAG = QMagenConf.Setting.PLOTFLAG;
%         catch
%             PLOTFLAG = 0;
%         end
%         
%         global EVOFLAG
%         EVOFLAG = 0;
%         
%         global SAVEFLAG
%         SAVEFLAG = QMagenConf.Setting.SAVEFLAG;
%         
%         global SAVENAME
%         SAVENAME = [QMagenConf.Setting.SAVENAME, '_', QMagenConf.Config.TStr, '/'];
%         
%         global SAVE_COUNT
%         if isempty(SAVE_COUNT)
%             SAVE_COUNT = 1;
%         end
%         
%         global MIN_LOSS_VAL
%         if isempty(MIN_LOSS_VAL)
%             MIN_LOSS_VAL = Inf;
%         end
%     case {'CALC-Cm', 'CALC-Chi'}
% end

% //Common parameters for many-body solvers: iLTRG, XTRG, ED, etc.
Para.Model = QMagenConf.ModelParaValue;
Para.Field = QMagenConf.Field;
Para.UnitCon = GetUnitCon(Para);
Para.beta_max = 1 / (K_min / Para.UnitCon.T_con);
[ Para ] = GetField( Para );

% //Import runtime parameters for each solver
Para = ImportMBSolverPara( Para );

switch QMagenConf.ModelConf.gFactor_Type
    case 'xyz'
        Para.g_fec.g = QMagenConf.ModelParaValue.g;
        Para.g_fec.dir = QMagenConf.ModelConf.gFactor_Vec;
        checka = sum(Para.g_fec.dir{1} .* Para.g_fec.dir{2});
        checkb = sum(Para.g_fec.dir{1} .* Para.g_fec.dir{3});
        checkc = sum(Para.g_fec.dir{2} .* Para.g_fec.dir{3});
        if abs(checka * checkb * checkc) * (norm(Para.g_fec.dir{1}) - 1) * (norm(Para.g_fec.dir{2}) - 1) * (norm(Para.g_fec.dir{3}) - 1)> 1e-15
            error('The directions of g factor are not perpendicular to each other or are not normalized!');
        end
    case 'dir'
        Para.g_fec.g = [QMagenConf.ModelParaValue.g, QMagenConf.ModelParaValue.g, QMagenConf.ModelParaValue.g];
        Para.g_fec.dir{1} = [1,0,0];
        Para.g_fec.dir{2} = [0,1,0];
        Para.g_fec.dir{3} = [0,0,1];
end
end