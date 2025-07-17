function FileName = GetFileName(Para)
FileName = ['Model=', Para.IntrcMap_Name(10:end), '-'];

fields = fieldnames(Para.Model);
for it = 1:length(fields)
    FileName = [FileName, fields{it}, '=', mat2str(eval(['Para.Model.', fields{it}])), '-'];
end

FileName = [FileName, 'h', '=', mat2str(Para.Field.h, 2), '-'];

fields = fieldnames(Para.Geo);
for it = 1:length(fields)
    FileName = [FileName, fields{it}, '=', num2str(eval(['Para.Geo.', fields{it}])), '-'];
end

if isfield(Para, 'saveInfo')
    FileName = [FileName, 'D=', num2str(Para.MCrit), '-', Para.saveInfo, '.mat'];
else
FileName = [FileName, 'D=', num2str(Para.MCrit), '.mat'];
end

end