function [ loss, Rslt ] = loss_func_chi( QMagenConf, Trange,  Chi_data, loss_type )
loss = 0;
Chi_data = sortrows(Chi_data, -1);
T_exp = Chi_data(:,1);
chi_exp = Chi_data(:,2);

T_min = max(min(T_exp), Trange(1));
T_max = min(max(T_exp), Trange(2));

ThDQ = 'Chi';

[~, Rslt] = GetResult(QMagenConf, 0.9 * T_min, ThDQ);
T = Rslt.T_l;
chi = Rslt.Chi_l;

switch QMagenConf.LossConf.IntSet
    case 'Int2Exp'
        for i = 1:length(T)
            if T(end) < T_min
                T(end) = [];
                chi(end) = [];
            else
                break;
            end
        end
        
        while 1
            if T(1) > T_max
                T(1) = [];
                chi(1) = [];
            else
                break;
            end
        end
        
        chi_int = interp1(T_exp, chi_exp, T);
        chi_f_exp = chi_int;
        chi_f_sim = chi;
        T_f = T;
        
    case 'Int2Sim'
        while 1
            if T_exp(1) > min(T_max, max(T))
                T_exp(1) = [];
                chi_exp(1) = [];
            else
                break;
            end
        end
        
        for i = 1:length(T)
            if T_exp(end) < T_min
                T_exp(end) = [];
                chi_exp(end) = [];
            else
                break;
            end
        end
        
        chi_int = interp1(T, chi, T_exp);
        chi_f_exp = chi_exp;
        chi_f_sim = chi_int;
        T_f = T_exp;
        
end

switch loss_type
    case 'abs-err'
        for i = 1:length(T_f)
            loss = loss + (chi_f_exp(i) - chi_f_sim(i))^2; 
        end
        loss = loss / max(chi_f_exp)^2;
    case 'rel-err'
        for i = 1:length(T)
            loss = loss + ((chi_f_exp(i) - chi_f_sim(i))/chi_f_sim(i))^2; 
        end
end

loss = loss/length(T);
end

