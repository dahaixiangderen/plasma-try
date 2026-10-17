function u = sps_unit(v)
%SPS_UNIT  单位化向量;当模长近似为 0 时返回零向量(不报错)
%
%   用于场方向、垂直基等场合,避免除零产生 NaN/Inf。

v = v(:);
n = norm(v);
if n <= eps
    u = zeros(size(v));
else
    u = v / n;
end

end
