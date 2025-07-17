function vec = GetPosi(Para, It_pos)

switch Para.IntrcMap_Name
    case {'IntrcMap_TLXXZ', 'IntrcMap_TLARX'}
        vec = GetTriPosi(Para, It_pos);
end
end

function vec = GetTriPosi(Para, It_pos)
Lx = Para.Geo.Lx;
Ly = Para.Geo.Ly;

a1 = [-1/2, -sqrt(3)/2];
a2 = [1, 0];
% here
for it = 1:length(It_pos)
    vec(it,:) = mod(It_pos(it)-1, Ly) * a1 + floor((It_pos(it)-1)/Ly) * a2; 
end
end
