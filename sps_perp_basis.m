function [e1, e2] = sps_perp_basis(kv, b0)
%SPS_PERP_BASIS  构造与 kv 垂直的右手正交基 (e1, e2, kv_hat)
%
%   kv : 参考方向(如波矢 k)。若为零向量则退化为 b0,再退化为 z 轴。
%   b0 : 备用参考方向(如背景磁场方向)
%
%   返回 e1, e2 均为单位向量且两两正交,e1 x e2 = kv_hat。

kh = sps_unit(kv);
if norm(kh) == 0
    kh = sps_unit(b0);
end
if norm(kh) == 0
    kh = [0; 0; 1];
end

tmp = [0; 0; 1];
if abs(dot(kh, tmp)) > 0.9
    tmp = [1; 0; 0];
end

e1 = sps_unit(cross(tmp, kh));
e2 = cross(kh, e1);

end
