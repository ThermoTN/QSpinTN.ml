function [ H ] = GetEDSMRslt( Para )

[ PMPO ] = Automata_SM( Para );

[ H ] = InitSM(PMPO, Para);
end

