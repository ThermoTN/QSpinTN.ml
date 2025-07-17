function savetdvpRslt(Para, Rslt)
FileName = GetFileName(Para);

save(['Rslt_DY/', FileName], 'Para', 'Rslt')
end