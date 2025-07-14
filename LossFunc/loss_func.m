function [ loss ] = loss_func(TStr, g, ParaVal, QMagenConf)

loss = 0;

TStrnow = datestr(now,'YYYYmmDD_HHMMSS');
len = length(TStrnow);

if len == length(TStr)
    QMagenConf.Config.TStr_log = TStrnow;
else
    QMagenConf.Config.TStr_log = [TStr(1:end-len), TStrnow];
end

ModelConf = QMagenConf.ModelConf;
CmData = QMagenConf.CmData;
ChiData = QMagenConf.ChiData;
LossConf = QMagenConf.LossConf;

wl = LossConf.WeightList;
loss_type = LossConf.Type;
loss_design = LossConf.Design;

if isempty(wl)
    wl = zeros(length(CmData) + length(ChiData), 1);
    wl(:) = 1;
end

ModelVal = struct();


for i = 1:1:length(ParaVal)
    ModelVal = setfield(ModelVal, QMagenConf.ModelConf.Para_Name{i}, ParaVal(i));
end
for i = 1:1:length(g)
    ModelVal = setfield(ModelVal, QMagenConf.ModelConf.gFactor_Name{i}, g(i));
end

QMagenConf.Config.fileID = GetFileID(QMagenConf, TStrnow, ModelVal);

% // Cm =======================================
RsltCv = cell(length(CmData), 1);
for i = 1:1:length(CmData)
    clear Field
    QMagenConf.Field.B = CmData(i).Info.Field;
    switch ModelConf.gFactor_Type
        case 'xyz'
            g_fec = g;
            
            if length(g_fec) ~= 3
                error('Not enough Lande g factor!');
            end
        case 'dir'
            g_fec = g(CmData(i).Info.g_info);
            
            check_flag = cross(QMagenConf.Field.B, ModelConf.gFactor_Vec{CmData(i).Info.g_info});
            if norm(QMagenConf.Field.B) ~= 0 && abs(norm(check_flag)) > 1e-10
                error('Effective Lande g factor is not along the field direction!');
            end
    end
    QMagenConf = GetModel(QMagenConf, g_fec, ParaVal);
    
    [lossp, RsltCv{i}] = loss_func_Cm(QMagenConf, CmData(i).Info.TRange, CmData(i).Data, loss_type);
    
    loss = loss + lossp * wl(i);
end

% // chi ======================================

RsltChi = cell(length(ChiData), 1);
for i = 1:1:length(ChiData)
    clear Field
    QMagenConf.Field.B = ChiData(i).Info.Field;
    if isnan(loss)
        break;
    end
    switch ModelConf.gFactor_Type
        case 'xyz'
            g_fec = g;
            
            if length(g_fec) ~= 3
                error('Not enough Lande g factor!');
            end
        case 'dir'
            g_fec = g(CmData(i).Info.g_info);
            
            check_flag = cross(QMagenConf.Field.B, ModelConf.gFactor_Vec{CmData(i).Info.g_info});
            if norm(QMagenConf.Field.B) ~= 0 && abs(norm(check_flag)) > 1e-10
                error('Effective Lande g factor is not along the field direction!');
            end
    end
    QMagenConf = GetModel(QMagenConf, g_fec, ParaVal);
    
    [lossp, RsltChi{i}] = loss_func_chi(QMagenConf, ChiData(i).Info.TRange, ChiData(i).Data, loss_type);
    
    loss = loss + lossp * wl(i + length(CmData));
end

switch loss_design
    case 'native'
    case 'log'
        loss = log10(loss);
    otherwise
        keyboard;
end

str = GetSaveFileName(QMagenConf, TStrnow, ModelVal);
if QMagenConf.Setting.SAVEFLAG == 1 && length(TStr) == length(TStrnow)
    save(str, 'ModelVal', 'RsltCv', 'RsltChi');
elseif QMagenConf.Setting.SAVEFLAG == 1 && length(TStr) ~= length(TStrnow)
    save(str, 'ModelVal', 'RsltCv', 'RsltChi');
end

end

function str = GetSaveFileName(QMagenConf, TStrnow, ModelVal)
str = ['../Tmp/', QMagenConf.Setting.SAVEFILENAME, QMagenConf.Config.TStr, '/'];
str = [str, 'Rslt'];
ModelName = fieldnames(ModelVal);
Value = struct2array(ModelVal);
for It = 1:length(ModelName)
    str = [str, '-', ModelName{It}, '=', num2str(Value(It))];
end
str = [str, '-time-', TStrnow, '.mat'];
end

function fileID = GetFileID(QMagenConf, TStrnow, ModelVal)
str = ['../Tmp/', QMagenConf.Setting.SAVEFILENAME, QMagenConf.Config.TStr, '/'];
str = [str, 'Log'];
ModelName = fieldnames(ModelVal);
Value = struct2array(ModelVal);
for It = 1:length(ModelName)
    str = [str, '-', ModelName{It}, '=', num2str(Value(It))];
end
str = [str, '-time-', TStrnow, '.log'];
fileID = fopen(str, 'a');
end