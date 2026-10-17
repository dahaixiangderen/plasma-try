function ok = sps_selftest()
%SPS_SELFTEST  物理与数值自检(可在无图形环境下运行)
%
%   ok = SPS_SELFTEST()   返回是否全部通过(逻辑值)
%
%   检验内容:
%     1) 均匀磁场纯回旋:动能、磁矩守恒,回旋周期/半径符合解析值
%     2) E x B 漂移:模拟得到的导心漂移等于 E/B
%     3) 偶极场:磁矩守恒、漂移方向(质子向西)、漂移速度与解析量级一致
%     4) 磁镜:损失锥内外粒子的捕获/逃逸行为
%     5) 各场模型的 div B ~ 0(数值检验)
%     6) 感应电场:E = -(omega/k) khat x dB 的量级检验
%     7) 相对论推进:纯磁场下 gamma 守恒
%     8) 全部预设至少能跑通若干步

fprintf('\n================= SPS_SELFTEST 开始 =================\n');
c = sps_const();
npass = 0; nfail = 0;

% ======================================================================
% 1) 均匀场纯回旋
% ======================================================================
cfg = sps_default_config();
cfg.verbose = false; cfg.plot = false;
cfg.species = 'proton';
cfg.B0 = [0;0;50e-9];
cfg.Emodel = 'zero';
cfg.init.speed = 1e5; cfg.init.alpha = 90; cfg.init.gyrophase = 0;
cfg.numerics.nper = 200; cfg.numerics.ngyro = 40;
out = sps_run(cfg);
s = out.summary;
rg_theory = out.cfg.scales.v0/out.cfg.scales.omega_c;
r_all = out.x - mean(out.x, 2);
rg_meas = mean(sqrt(sum(r_all.^2, 1)));
[npass, nfail] = chk('均匀场: |dKE|/KE < 1e-9', s.dKE_rel < 1e-9, s.dKE_rel, npass, nfail);
[npass, nfail] = chk('均匀场: |dmu|/mu < 1e-9', s.dmu_rel < 1e-9, s.dmu_rel, npass, nfail);
[npass, nfail] = chk('均匀场: 回旋周期误差 < 1e-3', ...
    abs(s.Tc_measured - s.Tc)/s.Tc < 1e-3, abs(s.Tc_measured-s.Tc)/s.Tc, npass, nfail);
[npass, nfail] = chk('均匀场: 回旋半径(均方根)误差 < 3%', ...
    abs(rg_meas - rg_theory)/rg_theory < 0.03, abs(rg_meas-rg_theory)/rg_theory, npass, nfail);

% ======================================================================
% 2) E x B 漂移
% ======================================================================
cfg = sps_default_config();
cfg.verbose = false; cfg.plot = false;
cfg.B0 = [0;0;50e-9];
cfg.Emodel = 'uniform'; cfg.E0 = [1e-3;0;0];
cfg.init.speed = 1e5; cfg.init.alpha = 90;
cfg.numerics.nper = 200; cfg.numerics.ngyro = 40;
out = sps_run(cfg);
vExB_theory = norm(cross(cfg.E0, cfg.B0))/norm(cfg.B0)^2;
vExB_meas   = norm(out.summary.vgc_measured);
[npass, nfail] = chk('E×B 漂移: 与 E/B 相对误差 < 2%', ...
    abs(vExB_meas - vExB_theory)/vExB_theory < 0.02, ...
    abs(vExB_meas-vExB_theory)/vExB_theory, npass, nfail);
% 漂移方向应为 E x B 方向
dir_theory = cross(cfg.E0, cfg.B0)/norm(cross(cfg.E0,cfg.B0));
dir_meas   = out.summary.vgc_measured/norm(out.summary.vgc_measured);
[npass, nfail] = chk('E×B 漂移: 方向正确(与 E×B 夹角 < 5 度)', ...
    acosd(min(max(dot(dir_theory,dir_meas),-1),1)) < 5, ...
    acosd(min(max(dot(dir_theory,dir_meas),-1),1)), npass, nfail);

% ======================================================================
% 3) 偶极场:磁矩守恒 + 漂移方向
% ======================================================================
cfg = sps_presets('dipole_belt');
cfg.verbose = false; cfg.plot = false;
cfg.numerics.ngyro = 600;
out = sps_run(cfg);
s = out.summary;
[npass, nfail] = chk('偶极场: 回旋平均 |dmu|/mu < 3e-2(绝热性)', ...
    s.dmu_smooth_rel < 3e-2, s.dmu_smooth_rel, npass, nfail);
fprintf('         [info] 瞬时 |dmu|/mu = %.4g ; 平滑 mu 均值 = %.6g J/T\n', ...
        s.dmu_rel, s.mu_smooth_mean);
[npass, nfail] = chk('偶极场: 瞬时 |dmu|/mu < 0.25(含 O(r_g*gradB/B) 摆动)', ...
    s.dmu_rel < 0.25, s.dmu_rel, npass, nfail);
[npass, nfail] = chk('偶极场: 质子向西漂移(v_gc,y < 0)', s.vgc_measured(2) < 0, ...
    s.vgc_measured(2), npass, nfail);
% 赤道偶极场的梯度-B 漂移量级(仅梯度项,曲率项在 alpha=90 时为 0)
d = sps_drifts(out.t(1), out.x(:,1), out.v(:,1), out.cfg);
v_gB = norm(d.vgradB);
[npass, nfail] = chk('偶极场: 漂移速度与梯度-B解析量级相符(0.5~2 倍)', ...
    v_gB > 0 && norm(s.vgc_measured)/v_gB > 0.5 && norm(s.vgc_measured)/v_gB < 2, ...
    norm(s.vgc_measured)/v_gB, npass, nfail);

% ======================================================================
% 4) 磁镜:损失锥内逃逸 / 损失锥外捕获
% ======================================================================
cfg = sps_presets('mirror_trap');
cfg.verbose = false; cfg.plot = false;
cfg.numerics.ngyro = 400;

cfgT = cfg; cfgT.init.alpha = 45;         % 损失锥角 30 度 -> 应被捕获
outT = sps_run(cfgT);
[npass, nfail] = chk('磁镜: alpha=45 度被捕获(检测到反弹)', ...
    isfinite(outT.summary.bounce_period) && outT.summary.n_bounce >= 2, ...
    outT.summary.n_bounce, npass, nfail);

cfgL = cfg; cfgL.init.alpha = 10;         % 损失锥内 -> 应沿场线逃逸
outL = sps_run(cfgL);
[npass, nfail] = chk('磁镜: alpha=10 度逃逸(未检测到反弹)', ...
    ~isfinite(outL.summary.bounce_period) || outL.summary.n_bounce < 2, ...
    outL.summary.n_bounce, npass, nfail);
[npass, nfail] = chk('磁镜: 反弹周期为正且有限(捕获粒子)', ...
    outT.summary.bounce_period > 0, outT.summary.bounce_period, npass, nfail);

% ======================================================================
% 5) 各场模型的 div B ~ 0
%    每行: 模型, 采样点, 差分步长 h, 参考波数 kref, 模型参数
%    判据: |div B| < 1e-2 * kref * |B|  (即相对误差 < 1%)
% ======================================================================
tests = { ...
 'uniform',         [1e6;0;0],      1e3, 1e-6,      struct(); ...
 'linear_gradient', [1e6;0;0],      1e3, 1e-6,      struct('Bdir',[0;0;1],'Bmag',50e-9, ...
                                                       'grad_dir',[1;0;0],'grad_amp',0.5,'L',1e6); ...
 'dipole',          [4*c.Re;0;0],   1e5, 1/(4*c.Re), struct('moment',[0;0;-c.m_earth]); ...
 'harris',          [0;0;0.3e6],    1e4, 1e-6,      struct('Bx',20e-9,'Bn',3e-9,'Bg',1e-9,'L',1e6); ...
 'mirror',          [1e5;0;0.5e6],  1e4, 1e-6,      struct('Bmag',100e-9,'L',1e6,'ratio',4); ...
 'wave',            [0;0;0],        200, 2.5e-4,    struct('k',[0;0;2.5e-4],'omega',1e4, ...
                                                       'B1',5e-9,'polarization','RH'); ...
 'parker',          [c.AU;0;0],     1e8, 1/c.AU,    struct('r0',10*c.Rsun,'B_r0',1.386e-6,'Vsw',4e5); ...
 'turbulent',       [1e6;0;0],      80,  2*pi/1e4,  struct('Nmodes',16,'dB_rms',1e-9,'seed',3); ...
};

for i = 1:size(tests,1)
    cfg = sps_default_config();
    cfg.verbose = false; cfg.plot = false;
    cfg.Bmodel   = tests{i,1};
    cfg.Bparams  = tests{i,5};
    cfg.init.mode = 'cart';
    cfg.init.x0   = tests{i,2}(:);
    cfg.init.v0   = [0;0;0];
    cfg.numerics.ngyro = 2;          % 本项只检验场,不做长时间积分
    cfg.numerics.nper  = 200;
    cfg = sps_prepare(cfg);

    x0 = cfg.init.x0; h = tests{i,3}; kref = tests{i,4};
    divB = 0;
    for k = 1:3
        dx = zeros(3,1); dx(k) = h;
        [~, Bp1] = sps_fields(0, x0 + dx, cfg);
        [~, Bp2] = sps_fields(0, x0 - dx, cfg);
        divB = divB + (Bp1(k) - Bp2(k))/(2*h);
    end
    [~, Bc] = sps_fields(0, x0, cfg);
    ref = kref*norm(Bc);
    [npass, nfail] = chk(sprintf('div B ~ 0: %-16s', tests{i,1}), ...
        abs(divB) < 1e-2*ref, abs(divB)/max(ref,realmin), npass, nfail);
end

% ======================================================================
% 6) 感应电场量级
% ======================================================================
cfg = sps_default_config();
cfg.verbose = false; cfg.plot = false;
kk = 2.5e-4; om = 1e4; B1 = 5e-9;
cfg.Bmodel = 'wave';
cfg.B0 = [0;0;200e-9];
cfg.Bparams = struct('k',[0;0;kk],'omega',om,'B1',B1,'polarization','RH');
cfg.Emodel = 'inductive_wave';
[E, ~] = sps_fields(0.3, [1e5;2e5;3e4], cfg);
E_theory = (om/kk)*B1;
[npass, nfail] = chk('感应电场: |E| = (omega/k)*B1 误差 < 1%', ...
    abs(norm(E) - E_theory)/E_theory < 0.01, abs(norm(E)-E_theory)/E_theory, npass, nfail);

% ======================================================================
% 6b) 波场/湍流场的 Faraday 定律一致性: curl E + dB/dt ~ 0
%     以及理想 Alfven 波无平行电场: E . b0 = 0
% ======================================================================
ft = { ...
 'wave',      [1e5;2e5;3e4], 200,  5e-6, ...
     struct('k',[0;0;2.5e-4],'omega',1e4,'B1',5e-9,'polarization','RH'), 200e-9; ...
 'turbulent', [0;0;0],       1.6e4, 0.33, ...
     struct('Nmodes',16,'dB_rms',1.5e-9,'seed',11), 5e-9; ...
};
for i = 1:size(ft,1)
    cfg = sps_default_config();
    cfg.verbose = false; cfg.plot = false;
    cfg.Bmodel  = ft{i,1};
    cfg.Bparams = ft{i,5};
    cfg.B0      = [0;0;ft{i,6}];
    if strcmp(ft{i,1},'wave')
        cfg.Emodel = 'inductive_wave';
    else
        cfg.Emodel = 'turbulence_inductive';
    end
    cfg.init.mode = 'cart';
    cfg.init.x0   = ft{i,2}(:);
    cfg.init.v0   = [0;0;0];
    cfg.numerics.ngyro = 2; cfg.numerics.nper = 200;
    cfg = sps_prepare(cfg);

    x0 = cfg.init.x0; h = ft{i,3}; ht = ft{i,4}; t0 = 0.3;
    dE = zeros(3,3);
    for j = 1:3
        dx = zeros(3,1); dx(j) = h;
        [Ep, ~] = sps_fields(t0, x0 + dx, cfg);
        [Em, ~] = sps_fields(t0, x0 - dx, cfg);
        dE(:,j) = (Ep - Em)/(2*h);
    end
    curlE = [dE(3,2)-dE(2,3); dE(1,3)-dE(3,1); dE(2,1)-dE(1,2)];
    [~, Bp] = sps_fields(t0 + ht, x0, cfg);
    [~, Bm] = sps_fields(t0 - ht, x0, cfg);
    dBdt = (Bp - Bm)/(2*ht);

    resid = norm(curlE + dBdt);
    scale = max([norm(curlE), norm(dBdt), realmin]);
    [npass, nfail] = chk(sprintf('Faraday 一致: %-12s', ft{i,1}), ...
        resid < 5e-2*scale, resid/scale, npass, nfail);

    [E0, B0v] = sps_fields(t0, x0, cfg);
    b0h = sps_unit(B0v);
    par  = abs(dot(E0, b0h))/max(norm(E0), realmin);
    [npass, nfail] = chk(sprintf('无平行电场 E.b0 ~ 0: %-8s', ft{i,1}), ...
        par < 1e-10, par, npass, nfail);
end

% ======================================================================
% 7) 相对论推进:纯磁场下 gamma 守恒
% ======================================================================
cfg = sps_default_config();
cfg.verbose = false; cfg.plot = false;
cfg.species = 'electron'; cfg.relativistic = true;
cfg.B0 = [0;0;200e-9];
cfg.Emodel = 'zero';
cfg.init.mode = 'cart';
cfg.init.x0 = [0;0;0];
cfg.init.v0 = [9e7; 4e7; 1e8];
cfg.numerics.nper = 200; cfg.numerics.ngyro = 50;
out = sps_run(cfg);
dg = max(abs(out.gamma - out.gamma(1)))/out.gamma(1);
[npass, nfail] = chk('相对论: 纯磁场下 gamma 守恒(<1e-9)', dg < 1e-9, dg, npass, nfail);

% ======================================================================
% 8) 全部预设能跑通(先用少量步数)
% ======================================================================
list = sps_preset_list();
for i = 1:numel(list)
    try
        cfg = sps_presets(list(i).name);
        cfg.verbose = false; cfg.plot = false;
        cfg.numerics.ngyro = 30;
        o = sps_run(cfg);
        okrun = numel(o.t) >= 5 && all(isfinite(o.x(:))) && all(isfinite(o.v(:)));
        [npass, nfail] = chk(sprintf('预设可运行: %-18s', list(i).name), okrun, ...
                             numel(o.t), npass, nfail);
    catch ME
        [npass, nfail] = chk(sprintf('预设可运行: %-18s', list(i).name), false, ...
                             ME.message, npass, nfail);
    end
end

% ======================================================================
% 9) 波-粒相互作用:回旋共振应使投掷角明显变化
% ======================================================================
wavs = {'whistler_electron','emic_proton'};
for i = 1:numel(wavs)
    cfg = sps_presets(wavs{i});
    cfg.verbose = false; cfg.plot = false;
    cfg.numerics.ngyro = 300;
    o = sps_run(cfg);
    [npass, nfail] = chk(sprintf('波-粒共振: %-18s 投掷角变化 > 0.5 度', wavs{i}), ...
        o.summary.pitch_spread > 0.5, o.summary.pitch_spread, npass, nfail);
end

% ======================================================================
fprintf('\n================= 自检结果: %d 通过, %d 失败 =================\n', npass, nfail);
ok = (nfail == 0);

end

% ----------------------------------------------------------------------
function [npass, nfail] = chk(name, cond, val, npass, nfail)
if cond
    fprintf('  [PASS] %-46s (值 = %s)\n', name, fmt(val));
    npass = npass + 1;
else
    fprintf(2, '  [FAIL] %-46s (值 = %s)\n', name, fmt(val));
    nfail = nfail + 1;
end
end

% ----------------------------------------------------------------------
function s = fmt(v)
if ischar(v)
    s = v;
elseif isempty(v)
    s = 'empty';
elseif isscalar(v) && isnumeric(v) && isfinite(v)
    s = sprintf('%.6g', v);
elseif isscalar(v) && isnumeric(v)
    s = 'NaN/Inf';
else
    s = sprintf('%dx%d %s', size(v,1), size(v,2), class(v));
end
end
