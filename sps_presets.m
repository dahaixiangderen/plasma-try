function [cfg, list] = sps_presets(name)
%SPS_PRESETS  空间等离子体环境预设算例(不同初始条件 / 不同场位形)
%
%   list = SPS_PRESETS()          返回所有预设的名称与说明
%   cfg  = SPS_PRESETS('name')    返回该预设的配置结构体(未 prepare)
%
%   所有预设都给出物理量级真实的数值(地球偶极场、1 AU 太阳风、磁尾电流片、
%   辐射带、哨声波/EMIC 波等),并把"时间变化"和"空间弯曲"体现在场模型里。
%   具体数值可用 sps_config_menu 交互修改,或在脚本里直接覆盖字段。
%
%   预设列表:
%     gyro_basic          均匀 B + 匀强 E:回旋运动与 E x B 漂移(教学)
%     gradient_drift      空间线性梯度场:梯度-B 漂移
%     dipole_belt         地球偶极场(弯曲场线)+ 1 MeV 质子:梯度+曲率漂移
%     convection          偶极场 + 时变对流电场:磁层对流与极化漂移
%     solar_wind          1 AU 太阳风:Parker 螺旋磁场 + 对流电场
%     tail_current_sheet  磁尾 Harris 电流片:重联位形中的混沌轨道
%     mirror_trap         磁镜位形:回弹运动与损失锥
%     whistler_electron   哨声波 - 相对论电子回旋共振
%     emic_proton         EMIC 波 - 质子回旋共振
%     turbulence          Alfven 湍流中的粒子轨道(随机场)

c = sps_const();

if nargin < 1 || isempty(name)
    cfg  = [];
    list = sps_preset_list();      % 等价于 list = sps_preset_list()
    return;
end

cfg = sps_default_config();
cfg.tag = name;

switch lower(strtrim(name))

    % ==================================================================
    case 'gyro_basic'
        % 教学算例:均匀磁场中的回旋运动 + 匀强电场引起的 E x B 漂移。
        % B = 50 nT (+z), E = 1 mV/m (+x) -> E x B 漂移 = E/B = 20 km/s (-y)
        cfg.tag      = '均匀B中的回旋与E×B漂移';
        cfg.species  = 'proton';
        cfg.Bmodel   = 'uniform';
        cfg.B0       = [0; 0; 50e-9];
        cfg.Emodel   = 'uniform';
        cfg.E0       = [1e-3; 0; 0];
        cfg.init.mode      = 'gc';
        cfg.init.x0        = [0;0;0];
        cfg.init.speed     = 1e5;        % 100 km/s
        cfg.init.alpha     = 90;         % 垂直投掷角
        cfg.init.gyrophase = 0;
        cfg.numerics.nper  = 120;
        cfg.numerics.ngyro = 20;

    % ==================================================================
    case 'gradient_drift'
        % |B| 沿 x 线性增强:B_z = 50 nT*(1 + 0.5*x/L),L = 1000 km。
        % 梯度-B 漂移 v = m*v_perp^2/(2*q*B^3) * (B x grad B) ~ 0.5 km/s (+y)
        cfg.tag      = '线性梯度磁场中的梯度-B漂移';
        cfg.species  = 'proton';
        cfg.Bmodel   = 'linear_gradient';
        cfg.B0       = [0; 0; 50e-9];
        cfg.Bparams  = struct('Bdir',[0;0;1], 'Bmag',50e-9, ...
                              'grad_dir',[1;0;0], 'grad_amp',0.5, 'L',1e6);
        cfg.Emodel   = 'zero';
        cfg.init.mode      = 'gc';
        cfg.init.x0        = [0;0;0];
        cfg.init.speed     = 1e5;
        cfg.init.alpha     = 90;
        cfg.init.gyrophase = 0;
        cfg.numerics.nper  = 120;
        cfg.numerics.ngyro = 200;        % 足够长以看清漂移

    % ==================================================================
    case 'dipole_belt'
        % 地球偶极场(真实弯曲场线),1 MeV 质子置于 L = 4 赤道面。
        % 弯曲场线 -> 梯度-B + 曲率漂移 -> 绕地球方位角漂移(辐射带漂移)
        Re = c.Re;
        cfg.tag      = '地球偶极场 L=4 的 1 MeV 质子';
        cfg.species  = 'proton';
        cfg.Bmodel   = 'dipole';
        cfg.B0       = [0; 0; c.B_eq_earth/4^3];   % 仅作参考量级
        cfg.Bparams  = struct('moment',[0;0;-c.m_earth], 'origin',[0;0;0]);
        cfg.Emodel   = 'zero';
        cfg.init.mode      = 'gc';
        cfg.init.x0        = [4*Re; 0; 0];
        cfg.init.speed     = sps_v_from_E(1e6, c.mp);   % 1 MeV 质子 ~ 1.38e7 m/s
        cfg.init.alpha     = 90;
        cfg.init.gyrophase = 0;
        cfg.numerics.nper  = 50;         % 偶极场漂移周期很长,用较少步/周期
        cfg.numerics.ngyro = 4000;       % ~ 540 s,可见明显方位角漂移

    % ==================================================================
    case 'convection'
        % 磁层对流:L = 6 偶极场 + 晨昏对流电场(0.5 mV/m)与 100 s ULF 扰动。
        % E x B 漂移把粒子向太阳方向输运,同时偶极场给出梯度/曲率漂移。
        Re = c.Re;
        cfg.tag      = '磁层对流(L=6, 偶极场+时变电场)';
        cfg.species  = 'proton';
        cfg.Bmodel   = 'dipole';
        cfg.B0       = [0; 0; c.B_eq_earth/6^3];
        cfg.Bparams  = struct('moment',[0;0;-c.m_earth]);
        cfg.Emodel   = 'time_harmonic';        % 时变电场!
        cfg.E0       = [0; 5e-4; 0];           % 晨昏电场 0.5 mV/m
        cfg.Eparams  = struct('E1',2e-4, 'Edir',[0;1;0], 'omega',2*pi/100);
        cfg.init.mode      = 'gc';
        cfg.init.x0        = [6*Re; 0; 0];
        cfg.init.speed     = sps_v_from_E(1e3, c.mp);   % 1 keV 质子
        cfg.init.alpha     = 90;
        cfg.init.gyrophase = 0;
        cfg.numerics.nper  = 60;
        cfg.numerics.ngyro = 600;

    % ==================================================================
    case 'solar_wind'
        % 1 AU 太阳风:Parker 螺旋磁场(B_r=3 nT, B_phi=-3.1 nT, 螺旋角 ~45 度)
        % + 对流电场 E = -U x B (U = 400 km/s) -> E x B 漂移 = 太阳风速度
        cfg.tag      = '1 AU 太阳风(Parker螺旋场+对流电场)';
        cfg.species  = 'proton';
        cfg.Bmodel   = 'parker';
        cfg.B0       = [3e-9; -3.06e-9; 0];    % 1 AU 处参考值
        cfg.Bparams  = struct('r0',10*c.Rsun, 'B_r0',1.386e-6, ...
                              'Vsw',4e5, 'OmegaSun',c.OmegaSun, 'axis',[0;0;1]);
        cfg.Emodel   = 'convective';
        cfg.Eparams  = struct('U',[4e5;0;0]);  % 径向太阳风
        cfg.init.mode      = 'gc';
        cfg.init.x0        = [c.AU; 0; 0];
        cfg.init.speed     = 3e4;              % 质子热速度 30 km/s
        cfg.init.alpha     = 45;
        cfg.init.gyrophase = 0;
        cfg.init.vdrift    = [4e5; 0; 0];      % 随太阳风整体流动(400 km/s)
        cfg.numerics.nper  = 100;
        cfg.numerics.ngyro = 200;

    % ==================================================================
    case 'tail_current_sheet'
        % 磁尾 Harris 电流片:B_x = 20 nT*tanh(z/L),L = 1000 km,
        % 法向场 B_n = 3 nT(重联位形)。粒子在片中做 Speiser 型混沌运动。
        cfg.tag      = '磁尾Harris电流片(重联位形)';
        cfg.species  = 'proton';
        cfg.Bmodel   = 'harris';
        cfg.B0       = [0; 0; 3e-9];
        cfg.Bparams  = struct('Bx',20e-9, 'Bn',3e-9, 'Bg',0, 'L',1e6, ...
                              'nhat',[0;0;1], 'tangent',[1;0;0]);
        cfg.Emodel   = 'zero';
        cfg.init.mode      = 'gc';
        cfg.init.x0        = [0;0;0];
        cfg.init.speed     = 5e5;              % 500 km/s
        cfg.init.alpha     = 60;
        cfg.init.gyrophase = 0;
        cfg.numerics.nper  = 100;
        cfg.numerics.ngyro = 200;

    % ==================================================================
    case 'mirror_trap'
        % 有限长度轴对称磁镜:B_z = 100 nT*[1+(R-1)sin^2(pi z/2L)],R = 4,L = 1000 km。
        % 镜比 R = 4 -> 损失锥角 asin(sqrt(1/R)) = 30 度:
        %   alpha < 30 度沿场线逃逸;alpha > 30 度被捕获并来回反弹。
        cfg.tag      = '磁镜捕获(辐射带反弹运动)';
        cfg.species  = 'proton';
        cfg.Bmodel   = 'mirror';
        cfg.B0       = [0; 0; 100e-9];
        cfg.Bparams  = struct('Bmag',100e-9, 'L',1e6, 'ratio',4, 'axis',[0;0;1]);
        cfg.Emodel   = 'zero';
        cfg.init.mode      = 'gc';
        cfg.init.x0        = [0;0;0];
        cfg.init.speed     = 2e5;              % 200 km/s (~200 eV)
        cfg.init.alpha     = 45;               % > 30 度 -> 被捕获;对比 10/25 度
        cfg.init.gyrophase = 45;
        cfg.numerics.nper  = 100;
        cfg.numerics.ngyro = 500;

    % ==================================================================
    case 'whistler_electron'
        % 哨声波 - 电子回旋共振(类似磁层合声波与 ~50 keV 电子的相互作用)
        % B0 = 200 nT, f = 2 kHz(< f_ce = 5.6 kHz), lambda ~ 25 km,
        % 回旋共振条件 omega - k*v_par = Omega_ce  ->  v_par ~ 9e7 m/s
        % 该电子为相对论速度(0.42c),故启用相对论 Boris 推进。
        n  = 5e6;                       % 电子密度 [m^-3]
        B0 = 200e-9;
        wpe = sqrt(n*c.e^2/(c.eps0*c.me));
        wce = c.e*B0/c.me;
        f   = 2000;
        om  = 2*pi*f;
        kk  = sqrt(om*wpe^2/(c.c_light^2*wce));     % 平行哨声波色散
        cfg.tag      = '哨声波-电子回旋共振(相对论电子)';
        cfg.species  = 'electron';
        cfg.relativistic = true;
        cfg.Bmodel   = 'wave';
        cfg.B0       = [0; 0; B0];
        cfg.Bparams  = struct('k',[0;0;-kk], 'omega',om, 'B1',10e-9, ...
                              'polarization','LH', 'e1',[1;0;0], 'phase0',0);
        cfg.Emodel   = 'inductive_wave';
        cfg.init.mode      = 'gc';
        cfg.init.x0        = [0;0;0];
        cfg.init.vpar      = 9.0e7;
        cfg.init.vperp     = 9.0e7;
        cfg.init.gyrophase = 0;
        cfg.numerics.nper  = 100;
        cfg.numerics.ngyro = 1000;

    % ==================================================================
    case 'emic_proton'
        % EMIC 波 - 质子回旋共振:B0 = 100 nT, f = 0.5 Hz(< f_cp = 1.53 Hz),
        % 左旋偏振,V_A = 690 km/s,k = 4.55e-6 rad/m(波长 ~ 1380 km)
        % 共振条件 omega - k*v_par = -Omega_cp -> v_par ~ 2.8e6 m/s
        n  = 1e7;
        B0 = 100e-9;
        mm = c.mp;
        VA = B0/sqrt(c.mu0*n*mm);
        f  = 0.5;
        om = 2*pi*f;
        kk = om/VA;
        cfg.tag      = 'EMIC波-质子回旋共振';
        cfg.species  = 'proton';
        cfg.Bmodel   = 'wave';
        cfg.B0       = [0; 0; B0];]);
        cfg.Emodel   = 'inductive_wave';
        cfg.init.mode      = 'gc';
        cfg.Bparams  = struct('k',[0;0;kk], 'omega',om, 'B1',5e-9, ...
                              'polarization','LH', 'e1',[1;0;0
        cfg.init.x0        = [0;0;0];
        cfg.init.speed     = 3.23e6;      % ~ 54 keV
        cfg.init.alpha     = 30;          % v_par = 2.8e6 m/s 恰好共振
        cfg.init.gyrophase = 0;
        cfg.numerics.nper  = 100;
        cfg.numerics.ngyro = 1000;

    % ==================================================================
    case 'turbulence'
        % 均匀背景场 + 多模 Alfven 湍流(不可压、随机相位、可复现)。
        % 波长 2e6 ~ 1e8 m(约 10 ~ 500 个回旋半径),最快模式周期 ~41 s
        % 长于回旋周期 13 s,粒子来得及响应,也不会白白消耗时间步。
        cfg.tag      = 'Alfvén湍流中的粒子轨道';
        cfg.species  = 'proton';
        cfg.Bmodel   = 'turbulent';
        cfg.B0       = [0; 0; 5e-9];       % 行星际磁场 ~5 nT
        cfg.Bparams  = struct('Nmodes',48, 'k_min',2*pi/1e8, 'k_max',2*pi/2e6, ...
                              'slope',5/3, 'dB_rms',1.5e-9, 'n0',5e6, 'seed',7);
        cfg.Emodel   = 'turbulence_inductive';
        cfg.init.mode      = 'gc';
        cfg.init.x0        = [0;0;0];
        cfg.init.speed     = 1e5;          % 100 km/s
        cfg.init.alpha     = 60;
        cfg.init.gyrophase = 0;
        cfg.numerics.nper  = 100;
        cfg.numerics.ngyro = 200;

    % ==================================================================
    otherwise
        error('sps_presets:unknown', ...
              '未知预设 "%s"。可用预设见 SPS_PRESETS() 的列表输出。', name);
end

end
