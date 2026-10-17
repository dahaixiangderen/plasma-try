function [modes, VA] = sps_make_turbulence(cfg)
%SPS_MAKE_TURBULENCE  生成多模 Alfven 湍流场(不可压、横向、随机相位)
%
%   [modes, VA] = SPS_MAKE_TURBULENCE(cfg)
%
%   每个模式 i:
%       dB_i(x,t) = dB(i) * [ e1*cos(k.x - omega*t + phase)
%                             + e2*sin(k.x - omega*t + phase) ]
%   其中 k 方向随机,e1,e2 垂直于 k,故 sum 场严格无散度;
%   幅度按幂律谱 dB(k) ~ k^(-slope/2) 分配,归一化到总均方根 dB_rms。
%
%   色散关系取剪切 Alfven 波: omega_i = VA * (k_i . b0)
%   (b0 为背景磁场方向;k.b0 = 0 的模式是准静态结构,omega = 0)。
%   这样配套的感应电场 E = -VA*(b0 x dB) 恰好满足 Faraday 定律
%   (逐模验证:nabla x E = -dB/dt),并且 E . b0 = 0,
%   即理想 MHD Alfven 波没有平行电场,不会出现虚假的平行加速。
%
%   cfg.Bparams 中可用参数:
%     Nmodes  (默认 48)       模式数
%     k_min   (默认 2pi/1e8)  最小波数 [rad/m]
%     k_max   (默认 2pi/2e6)  最大波数 [rad/m]
%     slope   (默认 5/3)      能谱指数
%     dB_rms  (默认 0.1*|B0|) 总扰动幅度 [T]
%     n0      (默认 5e6)      等离子体数密度 [m^-3],用于算 Alfven 速度
%     VA      (可选)          直接指定 Alfven 速度 [m/s]
%     seed    (默认 cfg.seed)
%
%   提示: 建议让最短波长 >= 几个回旋半径,并且使最快模式的周期 >= 回旋周期
%   (omega = k*VA <= Omega_c 附近)。周期远短于回旋周期的模式,粒子来不及响应,
%   只会白费时间步长(见 sps_prepare 中的 dt 判据)。

c = sps_const();
Bp = cfg.Bparams;

N     = round(sps_getp(Bp, 'Nmodes', 48));
kmin  = sps_getp(Bp, 'k_min', 2*pi/1e8);
kmax  = sps_getp(Bp, 'k_max', 2*pi/2e6);
slope = sps_getp(Bp, 'slope', 5/3);
B0    = cfg.B0(:);
B0m   = norm(B0);
if B0m <= 0, B0m = 1e-9; end
dB_rms = sps_getp(Bp, 'dB_rms', 0.1*B0m);
n0     = sps_getp(Bp, 'n0', 5e6);
VA     = sps_getp(Bp, 'VA', []);

if isempty(VA)
    switch lower(cfg.species)
        case {'electron','e','e-'}
            mm = c.me;
        otherwise
            mm = c.mp;
    end
    VA = B0m / sqrt(c.mu0 * n0 * mm);      % Alfven 速度 [m/s]
end

seed = sps_getp(Bp, 'seed', sps_getp(cfg,'seed',20240501));
if ~isempty(seed) && isnumeric(seed) && isfinite(seed)
    rng(seed, 'twister');
end

% ---- 波数抽样(对数均匀) ----
kk = exp(log(kmin) + rand(1,N)*(log(max(kmax,kmin*1.001)) - log(kmin)));
amp = kk.^(-slope/2);
amp = amp / sqrt(sum(amp.^2)) * dB_rms;    % 归一化使 sqrt(sum dB^2) = dB_rms

modes = struct('k',{}, 'omega',{}, 'dB',{}, 'e1',{}, 'e2',{}, 'phase',{}, 'kn',{});
b0h = sps_unit(B0);

for i = 1:N
    % 随机方向
    v = randn(3,1);
    while norm(v) < 1e-6
        v = randn(3,1);
    end
    kh = v/norm(v);
    [e1, e2] = sps_perp_basis(kh, [0;0;1]);
    modes(i).k     = kh*kk(i);
    modes(i).kn    = kk(i);
    modes(i).omega = kk(i)*VA*dot(kh, b0h);   % 剪切 Alfven 波色散 omega = VA*k_par
    modes(i).dB    = amp(i);
    modes(i).e1    = e1;
    modes(i).e2    = e2;
    modes(i).phase = 2*pi*rand();
end

end
