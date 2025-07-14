function [Rslt] = GettanTRGCorr(rho, Rslt, Para, Op)

[V, ~] = eig(reshape(Op.SmOp, [2,2]));
Sm = V(:,2) * V(:,1)';
Sp = Sm';

if length(Rslt.beta) == 2
    
    IntrMap = eval([Para.IntrcMap_Name, '(Para)']);
    for loc = 1:Para.L
        cnt = 1;
        clear IntrMPO
        minloc = Para.L;
        maxloc = 1;
        for it = 1:length(IntrMap)
            if ismember(loc, IntrMap(it).site)
                IntrMPO(cnt) = IntrMap(it);
                cnt = cnt + 1;
                minloc = min([minloc, IntrMap(it).site]);
                maxloc = max([maxloc, IntrMap(it).site]);
            end
        end
        
        for it = 1:length(IntrMPO)
            IntrMPO(it).site = IntrMPO(it).site - minloc + 1;
        end
        H = AutomataInit_MPO(IntrMPO, Para, (maxloc-minloc+1));
        H1 = H;
        H2 = H;
        if length(size(H1.A{loc - minloc+1})) == 3
            H1.A{loc - minloc+1} = contract(H1.A{loc - minloc+1}, 2, Sm, 2, [1,3,2]);
        else
            H1.A{loc - minloc+1} = contract(H1.A{loc - minloc+1}, 3, Sm, 2, [1,2,4,3]);
        end
        H2.A{loc - minloc+1} = contract(H2.A{loc - minloc+1}, length(size(H1.A{loc - minloc+1})), Sm, 1);
        Hcom.A = DircSumMPO(H1.A, H2.A, [1, -1]);
        try
            [Hcom.A, lgnorm] = CompressH(Hcom.A);
            Hcom.A{1} = exp(lgnorm) * Hcom.A{1};
            Hcom.lgnorm = 0;
        end
        Op_Hcom(loc).Hcom = Hcom;
        Op_Hcom(loc).minloc = minloc;
        Op_Hcom(loc).maxloc = maxloc;
    end
    str = num2str(randi(1000000));
    Rslt.Filestr = str;
    save(['CorrOpTmp/', str, '.mat'], 'Op_Hcom');
else
    load(['CorrOpTmp/', Rslt.Filestr, '.mat']);
end
EnVL = cell([Para.L-1, 1]);
for it = 1:Para.L-1
    if it == 1
        T = contract(rho.A{1}, [2,3], conj(rho.A{1}), [2,3]);
        EnVL{it} = T;
    else
        T = contract(T, 1, rho.A{it}, 1);
        T = contract(T, [1,3,4], conj(rho.A{it}), [1,3,4]);
        EnVL{it} = T;
    end
end

for it = 1:Para.L
    Tmploc = contract(rho.A{it}, length(size(rho.A{it})), conj(Sp), 2);
    minloc = Op_Hcom(it).minloc;
    for itj = Op_Hcom(it).minloc:Op_Hcom(it).maxloc
        if itj == Op_Hcom(it).minloc
            if itj == 1
                T = contract(rho.A{itj}, 2, Op_Hcom(it).Hcom.A{itj - minloc + 1}, 3);
            else
                T = contract(EnVL{itj-1}, 1, rho.A{itj}, 1);
                T = contract(T, 3, Op_Hcom(it).Hcom.A{itj - minloc + 1}, 3);
            end
            if itj == it
                if itj == 1
                    T = contract(T, [2,4], conj(Tmploc), [3,2]);
                else
                    T = contract(T, [1,3,4], conj(Tmploc), [1,4,3]);
                end
            else
                if itj == 1
                    T = contract(T, [2,4], conj(rho.A{itj}), [3,2]);
                else
                    T = contract(T, [1,3,5], conj(rho.A{itj}), [1,4,3]);
                end
            end
        elseif itj == Op_Hcom(it).maxloc
            if itj == Para.L
                T = contract(T, 1, rho.A{itj}, 1);
                T = contract(T, [1,3], Op_Hcom(it).Hcom.A{itj - minloc + 1}, [1,3]);
            else
                T = contract(T, 1, rho.A{itj}, 1);
                T = contract(T, [1,4], Op_Hcom(it).Hcom.A{itj - minloc + 1}, [1,3]);
            end
            if itj == it
                if itj == Para.L
                    T = contract(T, [1,3,2], conj(Tmploc), [1,2,3]);
                    ObsJ(it) = T;
                else
                    T = contract(T, [1,3,4,2], conj(Tmploc), [1,2,3,4]);
                    ObsJ(it) = T;
                end
            else
                if itj == Para.L
                    T = contract(T, [1,3,2], conj(rho.A{itj}), [1,2,3]);
                    ObsJ(it) = T;
                else
                    T = contract(T, [1,2,4,3], conj(rho.A{itj}), [1,2,3,4]);
                    ObsJ(it) = T;
                end
            end
        else
            T = contract(T, 1, rho.A{itj}, 1);
            T = contract(T, [1,4], Op_Hcom(it).Hcom.A{itj - minloc+1}, [1,4]);
            if itj == it
                T = contract(T, [1,5,3], conj(Tmploc), [1,3,4]);
            else
                T = contract(T, [1,5,3], conj(rho.A{itj}), [1,3,4]);
            end
        end
    end
end
ObsJ = real(ObsJ);
for it = 1:Para.L
    Tmp2 = contract(rho.A{it}, length(size(rho.A{it})), conj(Sp), 2);
    if it == 1
        Tmp1 = contract(rho.A{1}, 2, Sm, 2, [1,3,2]);
        ObsB(it) = contract(Tmp1, [1,2,3], conj(Tmp2), [1,2,3]);
    elseif it == Para.L
        Tmp1 = contract(rho.A{it}, 2, Sm, 2, [1,3,2]);
        Tmp1 = contract(EnVL{it-1}, 1, Tmp1, 1);
        ObsB(it) = contract(Tmp1, [1,2,3], conj(Tmp2), [1,2,3]);
    else
        Tmp1 = contract(rho.A{it}, 3, Sm, 2, [1,2,4,3]);
        Tmp1 = contract(EnVL{it-1}, 1, Tmp1, 1);
        ObsB(it) = contract(Tmp1, [1,2,3,4], conj(Tmp2), [1,2,3,4]);
    end
end
ObsB = real(ObsB) * norm(Para.Field.h);
if length(Rslt.beta) == 2
    Rslt.ObsJ{2} = ObsJ;
    Rslt.ObsB{2} = ObsB;
else
    Rslt.ObsJ{end + 1} = ObsJ;
    Rslt.ObsB{end + 1} = ObsB;
end
end