function [ Para ] = GetField( Para )

% try
%     Para.Field.h = Para.Field.h*1;
% catch
%     Para.Field.h = Para.Field.B ./ Para.UnitCon.h_con;
% end
Para.Field.h = norm(Para.Field.B) / Para.UnitCon.h_con;
Para.Field.hdir = Para.Field.B / norm(Para.Field.B);
end

