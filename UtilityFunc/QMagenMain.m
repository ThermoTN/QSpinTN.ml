function [ varargout ] = QMagenMain(QMagenConf, varargin)
% -------------------------------------------------
% Usage: 
%       callQSL(QMagenConf)
%       [Rslt_NU, Rslt_EXP] = callQSL(QMagenConf, 'Kmin', lowest set temperature)
% -------------------------------------------------
License_plot()
TStr = datestr(now,'YYYYmmDD_HHMMSS');
if ~exist('../Tmp','dir')
	mkdir('../Tmp');
end

switch QMagenConf.Config.Mode
    case 'OPT' 
        mkdir(['../Tmp/', QMagenConf.Setting.SAVEFILENAME, TStr]);
        
        QMagenConf.Config.TStr = TStr;
        disp(QMagenConf);
        OptFunc(TStr, QMagenConf)
        
    case 'ThDQ'
        mkdir(['../Tmp/', QMagenConf.Setting.SAVEFILENAME, TStr]);
        
        QMagenConf.Config.TStr = TStr;
        save(['../Tmp/', QMagenConf.Setting.SAVEFILENAME, TStr, '/configuration.mat'], 'QMagenConf');
        RsltCvNU = cell(length(varargin{2}), 1);
        RsltCv = cell(length(varargin{2}), 1);
        RsltChiNU = cell(length(varargin{2}), 1);
        RsltChi = cell(length(varargin{2}), 1);
        
        for loop_num = 1:1:length(varargin{2})
            ParaVal = zeros(1, length(QMagenConf.ModelConf.Para_Name));
            for i = 1:1:length(ParaVal)
                pos = find(strcmp(varargin, QMagenConf.ModelConf.Para_Name{i}), 1) + 1;
                try
                    ParaVal(i) = varargin{pos}(loop_num);
                catch
                    error('Not enough input arguemet! The value of %s is required!', ...
                            QMagenConf.ModelConf.Para_Name{i})
                end
            end
            
            g = zeros(1, QMagenConf.ModelConf.gFactor_Num);
            for i = 1:1:length(g)
                pos = find(strcmp(varargin, QMagenConf.ModelConf.gFactor_Name{i}), 1) + 1;
                try
                    g(i) = varargin{pos}(loop_num);
                catch
                    error('Not enough input arguemet! The value of %s is required!', ...
                            QMagenConf.ModelConf.gFactor_Name{i})
                end
            end
            [RsltCvNU{loop_num}, RsltCv{loop_num}, RsltChiNU{loop_num}, RsltChi{loop_num}] = GetThDQRslt( TStr, g, ParaVal );
        end
        varargout{1} = RsltCvNU;
        varargout{2} = RsltCv;
        varargout{3} = RsltChiNU;
        varargout{4} = RsltChi;
end