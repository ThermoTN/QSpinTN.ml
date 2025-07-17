function savetanTRGRslt(Para, Rslt)
FileName = GetFileName(Para);

save(['Rslt/', FileName], 'Para', 'Rslt')
end