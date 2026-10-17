function sps_demo(mode)
%SPS_DEMO  多算例演示:预设总览 + 参数扫描对比
%
%   SPS_DEMO()            全部演示
%   SPS_DEMO('presets')   只跑预设总览(结果汇总成表)
%   SPS_DEMO('sweep')     只跑参数扫描(投掷角/电场/弯曲场/时变场)
%
%   参数扫描是理解"不同初始条件 -> 不同轨道"的最快方式:
%     1) 磁镜:投掷角 10/25/45/70 度 -> 逃逸与捕获、反弹周期
%     2) 均匀场:电场 0/0.5/1/2 mV/m -> E x B 漂移速度
%     3) 弯曲 vs 直场:偶极场 vs 同强度均匀场 -> 曲率/梯度漂移
%     4) 场的时变快慢 omega/Omega_c = 0.05/0.5/3 -> 绝热性(磁矩是否守恒)

if nargin < 1 || isempty(mode)
    mode = 'all';
end
mode = lower(mode);

if any(strcmp(mode, {'all','presets'}))
    demo_presets();
end
if any(strcmp(mode, {'all','sweep'}))
    demo_sweep_alpha();
    demo_sweep_E();
    demo_sweep_curvature();
    demo_sweep_Bmod();
end

end

% ======================================================================
function demo_presets()
%DEMO_PRESETS  依次运行所有预设,汇总成一张对比表

fprintf('\n================= 预设算例总览 =================\n');
list = sps_preset_list();
res = struct('name',{}, 'tag',{}, 'B_nT',{}, 'Tc',{}, 'rg_km',{}, ...
             'mu',{}, 'vdrift',{}, 'vanal',{}, 'esc',{}, 'reason',{});
for i = 1:numel(list)
    cfg = sps_presets(list(i).name);
    cfg.verbose = false;
    cfg.plot = false;
    cfg.numerics.ngyro = min(cfg.numerics.ngyro, 300);   % 演示用,缩短时长
    fprintf('  [%2d/%2d] %-20s ...\n', i, numel(list), list(i).name);
    o = sps_run(cfg);
    s = o.summary;
    res(i).name   = list(i).name;
    res(i).tag    = cfg.tag;
    res(i).B_nT   = o.cfg.scales.B_now*1e9;
    res(i).Tc     = o.cfg.scales.Tc;
    res(i).rg_km  = o.cfg.scales.rg0/1e3;
    res(i).mu     = s.dmu_smooth_rel;
    res(i).vdrift = s.drift_speed_measured;
    res(i).vanal  = s.drift_speed_analytic;
    res(i).esc    = o.escaped;
    res(i).reason = o.esc_reason;
end

fprintf('\n%-20s %10s %10s %10s %12s %12s %12s\n', ...
        '预设', '|B|(nT)', 'T_c(s)', 'r_g(km)', '|dmu|/mu*', 'v_drift', 'v_解析');
fprintf('%s\n', repmat('-', 1, 94));
for i = 1:numel(res)
    fprintf('%-20s %10.3g %10.3g %10.4g %12.2e %12.4g %12.4g%s\n', ...
            res(i).name, res(i).B_nT, res(i).Tc, res(i).rg_km, ...
            res(i).mu, res(i).vdrift, res(i).vanal, ...
            ternary(res(i).esc, '  <提前终止>', ''));
end
fprintf('\n* |dmu|/mu 为回旋平均后的值(真正的绝热不变量)。\n');
fprintf('  v_drift = 模拟轨道端点给出的导心平均速度; v_解析 = 导心漂移解析合成。\n');
fprintf('  在薄电流片等 r_g 与场尺度可比的位形中,v_解析 的导心展开精度有限。\n');
fprintf('  捕获/反弹轨道(如 mirror_trap)的端点位移会因往返相消,该行对比不具意义。\n');
for i = 1:numel(res)
    if res(i).esc
        fprintf('  [提前终止] %-20s %s\n', res(i).name, res(i).reason);
    end
end

end

% ======================================================================
function demo_sweep_alpha()
%DEMO_SWEEP_ALPHA  磁镜中不同投掷角:逃逸 vs 捕获

alphas = [10 25 45 70];
fprintf('\n=== 扫描 1:磁镜中的投掷角 (损失锥角 30 度) ===\n');
figure('Name','扫描:磁镜投掷角','Color','w');
set(gcf,'Units','normalized','OuterPosition',[0.03 0.55 0.44 0.4]);
hold on; grid on;
leg = {};
for a = alphas
    cfg = sps_presets('mirror_trap');
    cfg.verbose = false; cfg.plot = false;
    cfg.init.alpha = a;
    cfg.numerics.ngyro = 400;
    o = sps_run(cfg);
    plot(o.t, o.pitch, 'LineWidth', 1.2);
    leg{end+1} = sprintf('\\alpha_0 = %d^o', a); %#ok<AGROW>
    if isfinite(o.summary.bounce_period)
        fprintf('  alpha=%2d 度: 被捕获,反弹周期 %.3g s, 投掷角范围 %.1f~%.1f\n', ...
                a, o.summary.bounce_period, o.summary.pitch_min, o.summary.pitch_max);
    else
        fprintf('  alpha=%2d 度: 落到损失锥内,沿场线逃逸\n', a);
    end
end
xlabel('时间 (s)'); ylabel('投掷角 (deg)'); ylim([0 180]);
title('磁镜:不同初始投掷角的轨道(损失锥 30^o)');
legend(leg, 'Location','best');

end

% ======================================================================
function demo_sweep_E()
%DEMO_SWEEP_E  均匀场中不同电场:E x B 漂移速度

Es = [0 0.5 1 2]*1e-3;
fprintf('\n=== 扫描 2:匀强电场大小 (B = 50 nT) ===\n');
figure('Name','扫描:电场与E×B漂移','Color','w');
set(gcf,'Units','normalized','OuterPosition',[0.5 0.55 0.44 0.4]);
hold on; grid on;
leg = {};
for E = Es
    cfg = sps_presets('gyro_basic');
    cfg.verbose = false; cfg.plot = false;
    cfg.E0 = [E; 0; 0];
    o = sps_run(cfg);
    plot(o.gc(1,:)/1e3, o.gc(2,:)/1e3, 'LineWidth', 1.2);
    leg{end+1} = sprintf('E = %.1f mV/m  (v_{E×B} = %.1f km/s)', ...
                         E*1e3, norm(o.summary.vgc_measured)/1e3); %#ok<AGROW>
    fprintf('  E = %.1f mV/m: 实测 E×B 漂移 = %.4g m/s (理论 %.4g m/s)\n', ...
            E*1e3, norm(o.summary.vgc_measured), E/norm(cfg.B0));
end
xlabel('x (km)'); ylabel('y (km)');
title('匀强电场下的导心轨迹:E×B 漂移');
legend(leg, 'Location','best'); axis equal;

end

% ======================================================================
function demo_sweep_curvature()
%DEMO_SWEEP_CURVATURE  弯曲场线(偶极场) vs 直线场(同 |B|)

fprintf('\n=== 扫描 3:弯曲场线 vs 直线场(相同 |B|) ===\n');
figure('Name','扫描:弯曲场与直场对比','Color','w');
set(gcf,'Units','normalized','OuterPosition',[0.03 0.06 0.44 0.4]);
hold on; grid on;

% (a) 偶极场:场线弯曲 -> 梯度+曲率漂移
cfg = sps_presets('dipole_belt');
cfg.verbose = false; cfg.plot = false;
cfg.numerics.ngyro = 1200;
o1 = sps_run(cfg);
plot((o1.gc(1,:)-o1.gc(1,1))/6.3712e6, (o1.gc(2,:)-o1.gc(2,1))/6.3712e6, 'LineWidth', 1.5);

% (b) 相同 |B| 的均匀直线场:alpha 相同,应只有 E×B(此处 E=0)
cfg = sps_default_config();
cfg.verbose = false; cfg.plot = false;
cfg.species = 'proton';
cfg.Bmodel  = 'uniform';
cfg.B0      = [0;0;3.12e-5/4^3];
cfg.init.mode = 'gc';
cfg.init.x0 = [0;0;0];
cfg.init.speed = sps_v_from_E(1e6, sps_const().mp);   % 同样 1 MeV 质子
cfg.init.alpha = 90;
cfg.numerics.ngyro = 1200;
o2 = sps_run(cfg);

% 画导心相对起点的位移,便于直接比较漂移
plot((o1.gc(1,:)-o1.gc(1,1))/6.3712e6, (o1.gc(2,:)-o1.gc(2,1))/6.3712e6, 'LineWidth', 1.5);
plot((o2.gc(1,:)-o2.gc(1,1))/6.3712e6, (o2.gc(2,:)-o2.gc(2,1))/6.3712e6, 'LineWidth', 1.5);

xlabel('x (R_E)'); ylabel('y (R_E)');
title('偶极场(弯曲)与均匀场(直线)中 1 MeV 质子的导心');
legend({sprintf('偶极场: |v_d| = %.3g m/s', norm(o1.summary.vgc_measured)), ...
        sprintf('均匀场: |v_d| = %.3g m/s', norm(o2.summary.vgc_measured))}, ...
       'Location','best');
axis equal;
fprintf('  偶极场(弯曲场线)导心漂移  = %.4g m/s\n', norm(o1.summary.vgc_measured));
fprintf('  均匀场(直线场线)导心漂移  = %.4g m/s (应接近 0)\n', norm(o2.summary.vgc_measured));

end

% ======================================================================
function demo_sweep_Bmod()
%DEMO_SWEEP_BMOD  场的时变快慢 -> 绝热性(磁矩守恒)是否成立
%
%   用**自洽**的平面波场:Bmodel='wave' 给出随时间的磁场扰动,
%   Emodel='inductive_wave' 按 Faraday 定律给出配套的感应电场。
%   扫描 omega/Omega_c:
%     omega << Omega_c : 场变化远慢于回旋 -> 绝热,磁矩守恒
%     omega ~  Omega_c : 场在一个回旋周期内显著变化 -> 非绝热,投掷角散射
%
%   (注: 若只想做"磁场强度随时间变化"的运动学测试,可用 cfg.Bmod;
%    它不含 Faraday 感应电场,因此不是自洽电磁场,磁矩会随 1/|B| 变化。)

c = sps_const();
ratios = [0.05 0.5 3];
fprintf('\n=== 扫描 4:场的时变快慢(omega/Omega_c)与绝热性 ===\n');
figure('Name','扫描:时变场与绝热性','Color','w');
set(gcf,'Units','normalized','OuterPosition',[0.5 0.06 0.44 0.4]);
hold on; grid on;
leg = {};
for r = ratios
    cfg = sps_default_config();
    cfg.verbose = false; cfg.plot = false;
    cfg.species = 'proton';
    B0  = 50e-9;
    Om  = c.e*B0/c.mp;              % Omega_c
    om  = r*Om;
    kk  = om/1e5;                   % 固定相速度 100 km/s(准静态极限)
    cfg.Bmodel  = 'wave';
    cfg.B0      = [0;0;B0];
    cfg.Bparams = struct('k',[0;0;kk], 'omega',om, 'B1',0.1*B0, ...
                         'polarization','RH', 'e1',[1;0;0]);
    cfg.Emodel  = 'inductive_wave';
    cfg.init.mode = 'gc';
    cfg.init.x0 = [0;0;0];
    cfg.init.speed = 1e5;           % 100 km/s
    cfg.init.alpha = 60;
    cfg.init.gyrophase = 0;
    cfg.numerics.nper  = 200;
    cfg.numerics.ngyro = 200;
    o = sps_run(cfg);
    plot(o.t/o.cfg.scales.Tc, o.pitch, 'LineWidth', 1.2);
    leg{end+1} = sprintf('\\omega/\\Omega_c = %-5.3g  |\\Delta\\mu|/\\mu = %.2e', ...
                         r, o.summary.dmu_smooth_rel); %#ok<AGROW>
    fprintf('  omega/Omega_c = %-5.3g : 投掷角 %.3f~%.3f 度, |dmu|/mu(回旋平均) = %.3e\n', ...
            r, o.summary.pitch_min, o.summary.pitch_max, o.summary.dmu_smooth_rel);
end
xlabel('t / T_c'); ylabel('投掷角 (deg)');
title('场时变快慢与绝热性(平面波 + 感应电场, \deltaB/B = 0.1)');
legend(leg, 'Location','best');

end

% ======================================================================
function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
