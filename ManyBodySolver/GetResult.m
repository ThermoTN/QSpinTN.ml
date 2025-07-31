function [Rslt, Rslt_exp_unit] = GetResult(QMagenConf, K_min, ThDQ)

Para = PassPara(QMagenConf, K_min);
Para.ThDQ = ThDQ;
Para.WorkingMode = QMagenConf.Config.Mode;
switch QMagenConf.Config.Mode
    case 'ThDQ'
        Para.fileID = 1;
        fprintf('\n')
        switch ThDQ
            case 'Cm'
                fprintf('Job started! \n%s\n', 'Heat capacity')
            case 'Chi'
                fprintf('Job started! \n%s\n', 'Magnetic susceptibility')
        end
        QMagenConf.Field.h = Para.Field.h;
        disp(QMagenConf);
        tstart = tic;
    case 'OPT'
        fileID = QMagenConf.Config.fileID;
        % fileID = fopen(['../Tmp/tmp_', Para.TStr, '/', Para.TStr_log, '.log'], 'a');
        Para.fileID = fileID;
        fprintf(fileID, '-------------------------------------------------------------------------------------------\n');
        switch ThDQ
            case 'Cm'
                fprintf(fileID, 'Job started! \n%s\n', 'Heat capacity');
            case 'Chi'
                fprintf(fileID, 'Job started! \n%s\n', 'Magnetic susceptibility');
        end
        fprintf(fileID, 'Parameters Value: \n');
        for i = 1:1:length(QMagenConf.ModelConf.Para_Name)
            if strcmp(QMagenConf.ModelConf.Para_Unit{i}, 'K')
                unit = 'Kelvin';
                CS = eval(['QMagenConf.ModelParaValue.', QMagenConf.ModelConf.Para_Name{i}]) * QMagenConf.ModelParaValue.ES;
            else
                unit = ['[', QMagenConf.ModelConf.Para_EnScale, ']'];
                CS = eval(['QMagenConf.ModelParaValue.', QMagenConf.ModelConf.Para_Name{i}]);
            end
            fprintf(fileID, '\t%-18s   |  %.3f %s\n', QMagenConf.ModelConf.Para_Name{i}, CS, unit);
        end
        fprintf(fileID, '\t%-18s   |  [%.3f, %.3f, %.3f] Tesla\n', 'Field', Para.Field.B);
        fprintf(fileID, '-------------------------------------------------------------------------------------------\n');
end

% //set up runtime parameters
% addpath('../ManyBodySolver/TensorFunction')
switch Para.ManyBodySolver
    case {'ED'}
%         addpath(genpath('../ManyBodySolver/ED'));
        Rslt = GetEDRslt(Para);
        
%         rmpath(genpath('../ManyBodySolver/ED'));
        
    case {'ED_SM'}
%         addpath(genpath('../ManyBodySolver/ED_SM'));
        Rslt = GetEDSMRslt(Para);
        
%         rmpath(genpath('../ManyBodySolver/ED_SM'));
    case {'ED_C'}
%         addpath(genpath('../ManyBodySolver/ED_C'));
        Rslt = GetEDCRslt(Para, ThDQ);
        
%         rmpath(genpath('../ManyBodySolver/ED_C'));
        
    case {'XTRG'}
%         addpath(genpath('../ManyBodySolver/XTRG'));
        
        Rslt = GetXTRGRslt(Para);
        
%         rmpath(genpath('../ManyBodySolver/XTRG'));
        
    case {'iLTRG'}
%         addpath(genpath('../ManyBodySolver/iLTRG'));
        
        Rslt = GetiLTRGRslt(Para, ThDQ);
        
%         rmpath(genpath('../ManyBodySolver/iLTRG'));
        
    case {'tanTRG'}
%         addpath(genpath('../ManyBodySolver/tanTRG'));
        
        Rslt = GettanTRGRslt(Para);
        Rslt.betaOri = Rslt.beta;
%         rmpath(genpath('../ManyBodySolver/tanTRG'));
        
    otherwise
        fprintf('Illegal mant-body solver! \n');
        keyboard;
end

Rslt.betaOri = Rslt.beta;
if strcmp(ThDQ, 'Cm') && (strcmp(Para.ManyBodySolver, 'XTRG') || strcmp(Para.ManyBodySolver, 'tanTRG'))
    [Rslt.CCHA, Rslt.betaCHA] = ThDQ_func( Rslt );
    Rslt.beta = Rslt.betaCHA;
    Rslt.Cm = Rslt.CCHA;
end

try
    Rslt_exp_unit.T_l = 1 ./ Rslt.beta * Para.UnitCon.T_con;
    Rslt_exp_unit.Cm_l = real(Rslt.Cm) * Para.UnitCon.Cm_con;
catch
    Rslt_exp_unit.Cm_l = zeros(size(Rslt.beta));
end

if norm(Para.Field.h) ~= 0 && strcmp(ThDQ, 'Chi')
    
    beta_list = 0.0025.*2.^(0:1:15).*2^0.1;
    for int = 0.3:0.2:0.9
        beta_list = [beta_list, 0.0025.*2.^(0:1:15).*2^int];
    end
    beta_list = sort(beta_list);
    loc1 = find(beta_list > Rslt.beta(1), 1);
    loc2 = find(beta_list > Rslt.beta(end), 1);
    if isempty(loc2)
        loc2 = length(beta_list);
    end
    beta = beta_list(loc1:loc2-1);
    Rslt.M = interp1(log(Rslt.betaOri), real(Rslt.M), log(beta), 'spline');
    Rslt.beta = beta;
    Rslt_exp_unit.T_l = 1 ./ Rslt.beta * Para.UnitCon.T_con;
    Rslt_exp_unit.Chi_l = real(Rslt.M) ./ Para.Field.h .* Para.UnitCon.Chi_con;
    Rslt_exp_unit.M = real(Rslt.M);
    
end


switch QMagenConf.Config.Mode
    case 'ThDQ'
        fprintf('Job finished! ')
        toc(tstart)
        fprintf('-------------------------------------------------------------------------------------------\n')
end

end