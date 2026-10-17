function d = sps_drifts(t, x, v, cfg)
%SPS_DRIFTS  计算导心漂移的各个解析分量(用于和模拟轨道对比)
%
%   d = SPS_DRIFTS(t, x, v, cfg)
%
%   采用导心理论的常用分解(非相对论形式):
%     平行运动    : v_par * b
%     电场漂移    : v_E   = E x B / B^2
%     梯度-B 漂移 : v_gB  = mu/(q B^3) * (B x grad|B|),  mu = m v_perp^2/(2B)
%     曲率漂移    : v_c   = m v_par^2/(q B^2) * (kappa x B),
%                   kappa = (b.grad)b 由中心差分求得
%     极化漂移    : v_p   = m/(q B^2) * dE_perp/dt
%
%   输出 d 的字段: vpar, vperp, mu, vExB, vgradB, vcurv, vparB, vpol, vtotal,
%                  Bmag, bhat, gradB, kappa, degenerate(是否接近零场)

x = x(:);
v = v(:);
q = cfg.q;
m = cfg.m;

[E, B, aux] = sps_fields(t, x, cfg);
Bn  = aux.Bmag;
bh  = aux.bhat;
h   = cfg.numerics.diff_h;

d = struct();
d.Bmag = Bn;
d.bhat = bh;
d.E = E;
d.B  = B;

if Bn < 1e-14 || norm(bh) == 0
    % 零场点:导心漂移公式不再适用
    d.degenerate = true;
    d.vpar = dot(v, [0;0;1]);
    d.vperp = norm(v - d.vpar*[0;0;1]);
    d.mu = NaN;
    d.vExB = [0;0;0]; d.vgradB = [0;0;0]; d.vcurv = [0;0;0];
    d.vparB = [0;0;0]; d.vpol = [0;0;0]; d.vtotal = v;
    d.gradB = [0;0;0]; d.kappa = [0;0;0];
    return;
end

% ---------- 场方向与梯度(中心差分) ----------
gradB = zeros(3,1);
kappa = zeros(3,1);
for i = 1:3
    dx = zeros(3,1);
    dx(i) = h;
    [~, Bp1] = sps_fields(t, x + dx, cfg);
    [~, Bp2] = sps_fields(t, x - dx, cfg);
    gradB(i) = (norm(Bp1) - norm(Bp2)) / (2*h);
    kappa   = kappa + bh(i) * (sps_unit(Bp1) - sps_unit(Bp2)) / (2*h);
end

% ---------- 速度分解 ----------
vpar     = dot(v, bh);
vperp_v  = v - vpar*bh;
vperp2   = sum(vperp_v.^2);
mu       = m*vperp2/(2*Bn);

% ---------- 各漂移分量 ----------
vExB   = cross(E, B) / Bn^2;
vgradB = (m*vperp2/(2*q*Bn^3)) * cross(B, gradB);
vcurv  = (m*vpar^2/(q*Bn^2))    * cross(kappa, B);
vparB  = vpar*bh;

% 极化漂移需要 dE/dt
ht = cfg.numerics.dt;
if isempty(ht) || ~(ht > 0)
    ht = cfg.scales.Tc/200;
end
[Ep, ~] = sps_fields(t + ht, x, cfg);
[Em, ~] = sps_fields(t - ht, x, cfg);
dEdt   = (Ep - Em)/(2*ht);
dEperp = dEdt - dot(dEdt, bh)*bh;
vpol   = (m/(q*Bn^2)) * dEperp;

d.degenerate = false;
d.vpar   = vpar;
d.vperp  = sqrt(vperp2);
d.mu     = mu;
d.vExB   = vExB;
d.vgradB = vgradB;
d.vcurv  = vcurv;
d.vparB  = vparB;
d.vpol   = vpol;
d.vtotal = vExB + vgradB + vcurv + vparB + vpol;
d.gradB  = gradB;
d.kappa  = kappa;

end
