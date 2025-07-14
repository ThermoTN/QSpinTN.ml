function [ loss, Rslt ] = loss_func_Cm( QMagenConf, Trange, C_data, loss_type )
loss = 0;
C_data = sortrows(C_data, -1);
T_exp = C_data(:,1);
C_exp = C_data(:,2);

T_min = max(min(T_exp), Trange(1));
T_max = min(max(T_exp), Trange(2));

ThDQ = 'Cm';

[~, Rslt] = GetResult(QMagenConf, 0.9 * T_min, ThDQ);
T = Rslt.T_l;
C = Rslt.Cm_l;

switch QMagenConf.LossConf.IntSet
    case 'Int2Exp'
        for i = 1:length(T)
            if T(end) < T_min
                T(end) = [];
                C(end) = [];
            else
                break;
            end
        end
        
        while 1
            if T(1) > T_max
                T(1) = [];
                C(1) = [];
            else
                break;
            end
        end
        
        C_int = interp1(T_exp, C_exp, T);
        C_f_exp = C_int;
        C_f_sim = C;
        T_f = T;
        
    case 'Int2Sim'
        while 1
            if T_exp(1) > min(T_max, max(T))
                T_exp(1) = [];
                C_exp(1) = [];
            else
                break;
            end
        end
        
        for i = 1:length(T)
            if T_exp(end) < T_min
                T_exp(end) = [];
                C_exp(end) = [];
            else
                break;
            end
        end
        
        C_int = interp1(T, C, T_exp);
        C_f_exp = C_exp;
        C_f_sim = C_int;
        T_f = T_exp;
        
end

switch loss_type
    case {'abs-err'}
        for i = 1:length(T_f)
            loss = loss + (C_f_exp(i) - C_f_sim(i))^2; 
        end
        loss = loss / max(C_f_exp)^2;
    case {'rel-err'}
        for i = 1:length(T_f)
            loss = loss + ((C_f_exp(i) - C_f_sim(i))/C_f_exp(i))^2; 
        end
end
loss = loss/length(T);
end

