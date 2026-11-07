function [CorSpin] = Correlate(Si, Sj, GS)
% Upper-triangular <Si_i Sj_j>. Both operators use (bra,ket) index order.
CorSpin = nan(length(GS), length(GS));
nGS = GS;
for it1 = 1:length(GS)-1
    if it1 == 1
        nGS{it1} = contract(nGS{it1}, 2, Si, 2);
    else
        nGS{it1} = contract(nGS{it1}, 3, Si, 2);
    end
    env = get_env(nGS, GS);
    nGS = GS;
    for it2 = it1+1:length(GS)
        if it2 == length(GS)
            loc_ten = contract(GS{it2}, 2, Sj, 2);
            tmp = contract(env{it2-1}, 1, loc_ten, 1);
            CorSpin(it1, it2) = contract(tmp, [1, 2], conj(GS{it2}), [1, 2]);
        else
            loc_ten = contract(GS{it2}, 3, Sj, 2);
            tmp = contract(env{it2-1}, 1, loc_ten, 1);
            tmp = contract(tmp, [1, 3], conj(GS{it2}), [1, 3]);
            CorSpin(it1, it2) = contract(tmp, [1, 2], eye(size(tmp, 1)), [1, 2]);
        end


    end
end
end

function env = get_env(nGS, GS)
env = cell(1, length(GS));

env{1} = contract(nGS{1}, 2, conj(GS{1}), 2);
for it = 2:length(GS)-1
    env{it} = contract(env{it-1}, 1, nGS{it}, 1);
    env{it} = contract(env{it}, [1, 3], conj(GS{it}), [1, 3]);
end

end
