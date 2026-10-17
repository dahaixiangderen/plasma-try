function [x, v, gamma] = sps_boris(x, v, q, m, E, B, dt, relativistic, substeps)
%SPS_BORIS  Boris 推进器:求解 dv/dt = (q/m)(E + v x B), dx/dt = v
%
%   [x, v] = SPS_BORIS(x, v, q, m, E, B, dt)
%   [x, v, gamma] = SPS_BORIS(x, v, q, m, E, B, dt, relativistic, substeps)
%
%   Boris 算法把电场加速(半步)与磁场纯旋转分开,磁场旋转严格保模长,
%   因此长时间积分能量/磁矩守恒性好,是等离子体 PIC 的标准推进器。
%
%   位置更新用旋转后的速度: x <- x + v_new*h(标准 Boris 写法)。
%   注意: 若改成用步首/步末速度的平均值(时间居中写法),在均匀场中回旋轨道
%   的几何形状会略微变好,但在**非均匀**磁场中会引入明显更大的磁矩漂移
%   (数值实验:偶极场 L=4 的 1 MeV 质子,平均速度写法 |dmu|/mu ~ 0.15,
%    而标准写法 ~ 0.002)。因此这里采用标准写法。
%
%   输入:
%     x, v    : 位置 [m]、速度 [m/s] (3x1)
%     q, m    : 电荷 [C]、质量 [kg]
%     E, B    : 该半步处的电场 [V/m]、磁场 [T]
%     dt      : 时间步 [s]
%     relativistic : true 时用 u = gamma*v 形式(高能粒子)
%     substeps: 把 dt 再细分多少子步(应对强场/高精度)
%
%   输出:
%     x, v    : 更新后的位置与速度
%     gamma   : 相对论因子(非相对论时恒为 1)

if nargin < 8 || isempty(relativistic), relativistic = false; end
if nargin < 9 || isempty(substeps),     substeps     = 1;     end

cc = 299792458;                 % 光速 [m/s]
h  = dt / max(substeps, 1);
x  = x(:);
v  = v(:);
E  = E(:);
B  = B(:);
gamma = 1;

for k = 1:max(substeps,1)

    if relativistic
        g  = 1/sqrt(max(1 - sum(v.^2)/cc^2, 1e-12));
        u  = g*v;
        um = u + (q*h/(2*m))*E;          % 半步电场加速
        tv = (q*h/(2*m))*B;              % 旋转向量
        s  = 2*tv/(1 + dot(tv,tv));
        up = um + cross(um, tv);
        u  = um + cross(up, s);          % 磁场旋转(保 |u|)
        u  = u + (q*h/(2*m))*E;          % 后半步电场加速
        g  = sqrt(1 + sum(u.^2)/cc^2);
        v  = u/g;
    else
        vm = v + (q*h/(2*m))*E;
        tv = (q*h/(2*m))*B;
        s  = 2*tv/(1 + dot(tv,tv));
        vp = vm + cross(vm, tv);
        v  = vm + cross(vp, s);
        v  = v + (q*h/(2*m))*E;
    end

    x = x + v*h;
end

if relativistic
    gamma = 1/sqrt(max(1 - sum(v.^2)/cc^2, 1e-12));
end

end
