function [Para] = InitKSpaceIntr(Para)


for it = 1:Para.L
    Intr(it).site = [it];
    Intr(it).operator = {Para.DY.QulocOp};
    Intr(it).J = exp(1i * Para.DY.k * it)/sqrt(Para.L);
end

Para.QuIntr = Intr;

end

