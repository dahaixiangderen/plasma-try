function cfg = sps_config_menu(cfg)
%SPS_CONFIG_MENU  交互式配置菜单(纯文本,所有数值都可手动输入)
%
%   cfg = SPS_CONFIG_MENU()      从"选择预设 / 默认配置"开始
%   cfg = SPS_CONFIG_MENU(cfg)   从给定配置开始修改
%
%   返回:修改后的 cfg;若用户选择取消则返回 []。
%
%   说明:数值输入框都接受表达式与向量,例如
%         2*pi/100 , 4*6.371e6 , [0 0 5e-9] , 1e5,0,0

if nargin < 1
    cfg = [];
end

% ---------------- 选择起点 ----------------
if isempty(cfg)
    list = sps_preset_list();
    fprintf('\n');
    fprintf('==================== 空间等离子体环境模拟:配置向导 ====================\n');
    fprintf('  0) 从默认配置开始(均匀 50 nT 磁场 + 100 km/s 质子)\n');
    for i = 1:numel(list)
        fprintf(' %2d) %-20s %-30s %s\n', i, list(i).name, list(i).title, list(i).physics);
    end
    fprintf('======================================================================\n');
    k = ask_num('选择起点编号', 0);
    if k >= 1 && k <= numel(list)
        cfg = sps_presets(list(k).name);
        fprintf('已载入预设 %s\n', list(k).name);
    else
        cfg = sps_default_config();
    end
end

% ---------------- 主菜单 ----------------
while true
    fprintf('\n================= 当前配置: %s =================\n', cfg.tag);
    fprintf('  粒子: %-16s  磁场: %-15s  电场: %-15s\n', ...
            cfg.species, cfg.Bmodel, cfg.Emodel);
    fprintf('  位置: [%.4g %.4g %.4g] m\n', cfg.init.x0(1), cfg.init.x0(2), cfg.init.x0(3));
    if strcmpi(cfg.init.mode,'gc')
        fprintf('  速率: %.4g m/s   投掷角: %.4g deg   (mode = gc)\n', ...
                cfg.init.speed, cfg.init.alpha);
    else
        fprintf('  速度: [%.4g %.4g %.4g] m/s  (mode = cart)\n', ...
                cfg.init.v0(1), cfg.init.v0(2), cfg.init.v0(3));
    end
    fprintf('--------------------------------------------------------------\n');
    fprintf('  1) 粒子种类 / 相对论开关\n');
    fprintf('  2) 磁场模型与空间位形(均匀/梯度/偶极/电流片/磁镜/波/Parker/湍流)\n');
    fprintf('  3) 磁场的时间变化(整体调制)\n');
    fprintf('  4) 电场模型与参数(零/均匀/时谐/圆偏振/对流/感应/平行/自定义)\n');
    fprintf('  5) 初始条件(位置、速率、投掷角、回旋相位、漂移)\n');
    fprintf('  6) 数值设置(步长/时长/步数/子步)\n');
    fprintf('  7) 输出与杂项(算例名/是否绘图/是否保存)\n');
    fprintf('  8) 打印配置摘要(试算一次,检查数值是否合理)\n');
    fprintf('  9) 运行模拟  <<<\n');
    fprintf('  0) 取消并返回\n');
    fprintf('--------------------------------------------------------------\n');

    sel = ask_num('请选择', 9);
    switch sel
        case 1, cfg = edit_species(cfg);
        case 2, cfg = edit_B(cfg);
        case 3, cfg = edit_Bmod(cfg);
        case 4, cfg = edit_E(cfg);
        case 5, cfg = edit_init(cfg);
        case 6, cfg = edit_numerics(cfg);
        case 7, cfg = edit_misc(cfg);
        case 8
            try
                sps_prepare(cfg);          % 只做校验并打印摘要
                fprintf('配置检查通过。\n');
            catch ME
                fprintf(2, '配置有问题: %s\n', ME.message);
            end
        case 9
            return;
        case 0
            cfg = [];
            return;
        otherwise
            fprintf('无效选项。\n');
    end
end

end

% ======================================================================
%                              粒子
% ======================================================================
function cfg = edit_species(cfg)
names = {'proton 质子 H+', 'electron 电子 e-', 'alpha He2+', 'O+ 氧离子', ...
         'custom 自定义(q,m 手动输入)'};
keys  = {'proton','electron','alpha','O+','custom'};
d = find(strcmpi(keys, cfg.species), 1);
if isempty(d), d = 1; end
i = ask_choice('粒子种类', names, d);
cfg.species = keys{i};
if strcmpi(cfg.species,'custom')
    cfg.q = ask_num('电荷 q (C,带符号,如 1.602e-19)', sps_getp(cfg,'q',1.602e-19));
    cfg.m = ask_num('质量 m (kg,如 1.673e-27)',      sps_getp(cfg,'m',1.673e-27));
end
cfg.relativistic = ask_yn('使用相对论 Boris 推进?', sps_getp(cfg,'relativistic',false));
end

% ======================================================================
%                              磁场
% ======================================================================
function cfg = edit_B(cfg)
names = { ...
 'uniform          均匀恒定磁场', ...
 'linear_gradient  |B| 沿某方向线性变化(梯度-B漂移)', ...
 'dipole           偶极场(弯曲场线,地球/行星)', ...
 'harris           Harris 电流片(磁尾重联)', ...
 'mirror           轴对称磁镜(捕获/损失锥)', ...
 'wave             平面波(哨声/EMIC,圆或线偏振)', ...
 'parker           Parker 螺旋(太阳风)', ...
 'turbulent        多模 Alfven 湍流', ...
 'custom           自定义函数 cfg.Bfun(t,x)'};
keys = {'uniform','linear_gradient','dipole','harris','mirror','wave','parker','turbulent','custom'};
d = find(strcmpi(keys, cfg.Bmodel), 1);
if isempty(d), d = 1; end
i = ask_choice('磁场模型', names, d);
cfg.Bmodel = keys{i};
Bp = cfg.Bparams;

switch cfg.Bmodel
    case 'uniform'
        cfg.B0 = ask_vec3('背景磁场 B0 (T, 例 50e-9=50nT)', cfg.B0);

    case 'linear_gradient'
        Bp.Bdir     = ask_vec3('磁场方向 Bdir', sps_getp(Bp,'Bdir',[0;0;1]));
        Bp.Bmag     = ask_num('参考点 |B| (T, 例 50e-9)', sps_getp(Bp,'Bmag',50e-9));
        Bp.grad_dir = ask_vec3('梯度方向 grad_dir', sps_getp(Bp,'grad_dir',[1;0;0]));
        Bp.grad_amp = ask_num('相对梯度强度 grad_amp(|B| 在 L 内变化倍数)', sps_getp(Bp,'grad_amp',0.5));
        Bp.L        = ask_num('梯度标长 L (m, 例 1e6=1000km)', sps_getp(Bp,'L',1e6));
        Bp.origin   = ask_vec3('场源原点 origin (m)', sps_getp(Bp,'origin',[0;0;0]));
        cfg.B0      = Bp.Bdir(:)*Bp.Bmag;

    case 'dipole'
        fprintf('提示: 地球偶极矩约 8.07e22 A*m^2(方向指向南,即 -z),\n');
        fprintf('      也可用 eq_field 给出赤道表面场(地球 3.12e-5 T)与 R_planet。\n');
        useq = ask_yn('用"赤道表面磁场+行星半径"给偶极矩?', ...
                      isfield(Bp,'eq_field') && ~isempty(sps_getp(Bp,'eq_field',[])));
        if useq
            Bp.eq_field = ask_num('赤道表面 |B| (T, 地球 3.12e-5)', sps_getp(Bp,'eq_field',3.12e-5));
            Bp.R_planet = ask_num('行星半径 R (m, 地球 6.3712e6)', sps_getp(Bp,'R_planet',6.3712e6));
            Bp.moment_dir = ask_vec3('偶极矩方向(地球取 [0 0 -1])', sps_getp(Bp,'moment_dir',[0;0;-1]));
            Bp = rmfield_safe(Bp, 'moment');
        else
            Bp.moment = ask_vec3('偶极矩 m (A*m^2, 例 [0 0 -8.07e22])', sps_getp(Bp,'moment',[0;0;-8.07e22]));
        end
        Bp.origin = ask_vec3('偶极中心 origin (m)', sps_getp(Bp,'origin',[0;0;0]));
        cfg.B0 = [0;0;3.12e-5/max(norm(cfg.init.x0)/6.3712e6,0.1)^3];   % 仅作参考

    case 'harris'
        Bp.L       = ask_num('电流片半厚度 L (m, 磁尾 ~1e6)', sps_getp(Bp,'L',1e6));
        Bp.Bx      = ask_num('反平行分量 Bx (T, 磁尾 ~2e-8)', sps_getp(Bp,'Bx',2e-8));
        Bp.Bn      = ask_num('法向分量 Bn (T, 重联 ~2e-9)', sps_getp(Bp,'Bn',2e-9));
        Bp.Bg      = ask_num('导向场 Bg (T)', sps_getp(Bp,'Bg',0));
        Bp.nhat    = ask_vec3('法向 nhat', sps_getp(Bp,'nhat',[0;0;1]));
        Bp.tangent = ask_vec3('切向 tangent', sps_getp(Bp,'tangent',[1;0;0]));
        Bp.origin  = ask_vec3('电流片中心 origin (m)', sps_getp(Bp,'origin',[0;0;0]));
        cfg.B0     = Bp.Bn(:)*sps_unit(Bp.nhat);

    case 'mirror'
        Bp.Bmag  = ask_num('中心 |B| (T, 例 100e-9)', sps_getp(Bp,'Bmag',100e-9));
        Bp.L     = ask_num('磁瓶半长 L (m, 场强在 |z|=L 处最大)', sps_getp(Bp,'L',1e6));
        Bp.ratio = ask_num('镜比 R = B_max/B_min (必须 >= 1, 例 4)', sps_getp(Bp,'ratio',4));
        Bp.axis  = ask_vec3('镜轴 axis', sps_getp(Bp,'axis',[0;0;1]));
        Bp.origin= ask_vec3('中心 origin (m)', sps_getp(Bp,'origin',[0;0;0]));
        cfg.B0   = Bp.axis(:)*Bp.Bmag;
        fprintf('  -> 损失锥角 alpha_loss = asin(sqrt(1/R)) = %.2f 度\n', ...
                asind(sqrt(1/max(Bp.ratio,1))));

    case 'wave'
        Bp.k   = ask_vec3('波矢 k (rad/m, 例 [0 0 2.5e-4])', sps_getp(Bp,'k',[0;0;2.5e-4]));
        Bp.omega = ask_num('角频率 omega (rad/s, 例 2*pi*2000)', sps_getp(Bp,'omega',2*pi*2000));
        Bp.B1  = ask_num('波振幅 B1 (T, 例 10e-9)', sps_getp(Bp,'B1',10e-9));
        polnames = {'RH','LH','linear'};
        dpol = find(strcmpi(polnames, sps_getp(Bp,'polarization','RH')), 1);
        if isempty(dpol), dpol = 1; end
        pol = ask_choice('偏振(相对波矢 k 的旋转方向)', ...
                         {'RH 右旋', 'LH 左旋', 'linear 线偏振'}, dpol);
        Bp.polarization = polnames{pol};
        Bp.phase0 = ask_num('初相位 phase0 (rad)', sps_getp(Bp,'phase0',0));
        Bp.e1  = ask_vec3('垂直基 e1(可留默认)', sps_getp(Bp,'e1',[]));
        if all(Bp.e1 == 0), Bp.e1 = []; end
        cfg.B0 = ask_vec3('背景磁场 B0 (T, 波场叠加在其上)', cfg.B0);

    case 'parker'
        Bp.r0       = ask_num('参考半径 r0 (m, 源面 ~10 R_sun=6.96e9)', sps_getp(Bp,'r0',6.96e9));
        Bp.B_r0     = ask_num('r0 处径向场 (T, 使 1AU 处 ~3nT 约 1.4e-6)', sps_getp(Bp,'B_r0',1.386e-6));
        Bp.Vsw      = ask_num('太阳风速度 (m/s, 400km/s=4e5)', sps_getp(Bp,'Vsw',4e5));
        Bp.OmegaSun = ask_num('太阳自转 (rad/s)', sps_getp(Bp,'OmegaSun',2.865e-6));
        Bp.axis     = ask_vec3('自转轴 axis', sps_getp(Bp,'axis',[0;0;1]));
        cfg.B0 = [3e-9;-3.06e-9;0];

    case 'turbulent'
        Bp.Nmodes = ask_num('模式数 Nmodes', sps_getp(Bp,'Nmodes',48));
        Bp.k_min  = ask_num('最小波数 k_min (rad/m, 长波)', sps_getp(Bp,'k_min',2*pi/1e8));
        Bp.k_max  = ask_num('最大波数 k_max (rad/m, 短波)', sps_getp(Bp,'k_max',2*pi/2e6));
        Bp.slope  = ask_num('能谱指数 slope(5/3=Kolmogorov)', sps_getp(Bp,'slope',5/3));
        Bp.dB_rms = ask_num('总扰动幅度 dB_rms (T)', sps_getp(Bp,'dB_rms',1.5e-9));
        Bp.n0     = ask_num('数密度 n0 (m^-3, 用于算 V_A)', sps_getp(Bp,'n0',5e6));
        Bp.seed   = ask_num('随机种子 seed', sps_getp(Bp,'seed',7));
        Bp        = rmfield_safe(Bp, 'modes');     % 参数变了要重新生成
        cfg.B0    = ask_vec3('背景磁场 B0 (T)', cfg.B0);

    case 'custom'
        fprintf('请在脚本中设置 cfg.Bfun = @(t,x) [Bx;By;Bz](单位 T)。\n');
        cfg.B0 = ask_vec3('背景磁场参考值 B0 (T)', cfg.B0);
end
cfg.Bparams = Bp;
end

% ======================================================================
function cfg = edit_Bmod(cfg)
fprintf('\n磁场整体时间调制: B(t) = B_model(t,x) * (1 + amp*sin(omega*t + phase))\n');
cfg.Bmod.amp   = ask_num('调制幅度 amp(0 = 无调制)', sps_getp(cfg.Bmod,'amp',0));
cfg.Bmod.omega = ask_num('角频率 omega (rad/s, 例 2*pi/600)', sps_getp(cfg.Bmod,'omega',0));
cfg.Bmod.phase = ask_num('初相位 phase (rad)', sps_getp(cfg.Bmod,'phase',0));
end

% ======================================================================
function cfg = edit_E(cfg)
names = { ...
 'zero                  E = 0', ...
 'uniform               均匀静电场 E0', ...
 'time_harmonic         E0 + E1*cos(omega t + phase)', ...
 'circular              垂直 B 的圆偏振电场(回旋加热)', ...
 'convective            对流电场 E = -U x B(太阳风/磁层对流)', ...
 'inductive_wave        波感应电场(Faraday,配合波场)', ...
 'parallel_wave         沿 B 的平行电场(朗道/平行加速)', ...
 'turbulence_inductive  湍流模感应电场', ...
 'custom                自定义函数 cfg.Efun(t,x)'};
keys = {'zero','uniform','time_harmonic','circular','convective', ...
        'inductive_wave','parallel_wave','turbulence_inductive','custom'};
d = find(strcmpi(keys, cfg.Emodel), 1);
if isempty(d), d = 1; end
i = ask_choice('电场模型', names, d);
cfg.Emodel = keys{i};
Ep = cfg.Eparams;

switch cfg.Emodel
    case 'zero'
        cfg.E0 = [0;0;0];
    case 'uniform'
        cfg.E0 = ask_vec3('静电场 E0 (V/m, 1mV/m = 1e-3)', cfg.E0);
    case 'time_harmonic'
        cfg.E0      = ask_vec3('恒定分量 E0 (V/m)', cfg.E0);
        Ep.E1       = ask_num('振荡幅度 E1 (V/m)', sps_getp(Ep,'E1',1e-4));
        Ep.Edir     = ask_vec3('振荡方向 Edir', sps_getp(Ep,'Edir',[1;0;0]));
        Ep.omega    = ask_num('角频率 omega (rad/s)', sps_getp(Ep,'omega',2*pi/100));
        Ep.phase    = ask_num('初相位 phase (rad)', sps_getp(Ep,'phase',0));
    case 'circular'
        cfg.E0   = ask_vec3('恒定分量 E0 (V/m)', cfg.E0);
        Ep.E1    = ask_num('旋转场幅度 E1 (V/m)', sps_getp(Ep,'E1',1e-4));
        Ep.omega = ask_num('角频率 omega (rad/s, 取 Omega_c 附近为回旋共振)', ...
                           sps_getp(Ep,'omega',1));
        Ep.phase = ask_num('初相位 phase (rad)', sps_getp(Ep,'phase',0));
        Ep.e1    = ask_vec3('垂直基 e1(可留默认)', sps_getp(Ep,'e1',[]));
        if all(Ep.e1 == 0), Ep.e1 = []; end
    case 'convective'
        Ep.U      = ask_vec3('等离子体流速 U (m/s, 太阳风 400km/s=4e5)', sps_getp(Ep,'U',[4e5;0;0]));
        Ep.Uamp   = ask_num('流速调制幅度(相对, 0=无)', sps_getp(Ep,'Uamp',0));
        Ep.Uomega = ask_num('流速调制频率 (rad/s)', sps_getp(Ep,'Uomega',0));
        fprintf('  -> E = -U x B,|E| = |U_perp|*|B|\n');
    case 'inductive_wave'
        fprintf('  -> 需要 Bmodel = wave 或 turbulent;E = -(omega/|k|) khat x dB\n');
    case 'parallel_wave'
        Ep.E1   = ask_num('平行电场幅度 E1 (V/m)', sps_getp(Ep,'E1',1e-5));
        Ep.phase = ask_num('初相位 phase (rad)', sps_getp(Ep,'phase',0));
    case 'turbulence_inductive'
        fprintf('  -> 需要 Bmodel = turbulent\n');
    case 'custom'
        fprintf('请在脚本中设置 cfg.Efun = @(t,x) [Ex;Ey;Ez](单位 V/m)。\n');
end
cfg.Eparams = Ep;
end

% ======================================================================
function cfg = edit_init(cfg)
fprintf('\n初始条件:\n');
modenames = {'gc','cart'};
dmode = find(strcmpi(modenames, cfg.init.mode), 1);
if isempty(dmode), dmode = 1; end
mode = ask_choice('初速度给定方式', ...
                  {'gc   用 速率/投掷角/回旋相位(推荐)', ...
                   'cart 直接给三维速度矢量'}, dmode);
cfg.init.mode = modenames{mode};

cfg.init.x0 = ask_vec3('初始位置 x0 (m)', cfg.init.x0);

if strcmpi(cfg.init.mode,'cart')
    cfg.init.v0 = ask_vec3('初始速度 v0 (m/s)', cfg.init.v0);
else
    fprintf('  说明: alpha 是与局地磁场的夹角;alpha=90 度表示纯垂直回旋。\n');
    cfg.init.speed     = ask_num('速率 |v| (m/s, 例 1e5=100km/s)', cfg.init.speed);
    cfg.init.alpha     = ask_num('投掷角 alpha (deg)', cfg.init.alpha);
    cfg.init.gyrophase = ask_num('回旋相位 (deg)', cfg.init.gyrophase);
    usevp = ask_yn('改用 v_par / v_perp 直接指定?', false);
    if usevp
        cfg.init.vpar  = ask_num('v_par (m/s)', 0);
        cfg.init.vperp = ask_num('v_perp (m/s)', 1e5);
    else
        cfg.init.vpar  = NaN;
        cfg.init.vperp = NaN;
    end
    cfg.init.vdrift = ask_vec3('附加漂移速度 vdrift (m/s)', cfg.init.vdrift);
end

% 位置的直观换算
c = sps_const();
r = norm(cfg.init.x0);
fprintf('  |x0| = %.4g m = %.4g km = %.4g R_E = %.4g AU\n', ...
        r, r/1e3, r/c.Re, r/c.AU);
end

% ======================================================================
function cfg = edit_numerics(cfg)
fprintf('\n数值设置(留空即保持默认):\n');
nu = cfg.numerics;
nu.nper = ask_num('每回旋周期步数 nper(推荐 50~200)', nu.nper);
nu.ngyro = ask_num('模拟时长 = 多少个初始回旋周期 ngyro', nu.ngyro);
dt = ask_num('直接指定时间步 dt (s, 0 = 按 nper 自动)', sps_getp(nu,'dt',0));
if dt > 0, nu.dt = dt; else, nu.dt = []; end
tm = ask_num('直接指定总时长 tmax (s, 0 = 按 ngyro 自动)', sps_getp(nu,'tmax',0));
if tm > 0, nu.tmax = tm; else, nu.tmax = []; end
nu.substeps     = ask_num('每步内 Boris 子步数(下限)', nu.substeps);
nu.adaptive     = ask_yn('局地磁场明显变强时自动增加子步(安全网)?', sps_getp(nu,'adaptive',true));
nu.substeps_trigger = ask_num('触发细分的倍数阈值(局地场/参考场)', sps_getp(nu,'substeps_trigger',2));
nu.substeps_max = ask_num('自适应子步上限', sps_getp(nu,'substeps_max',200));
nu.store_every  = ask_num('存储间隔 store_every(省内存)', nu.store_every);
cfg.numerics = nu;
end

% ======================================================================
function cfg = edit_misc(cfg)
cfg.tag      = ask_str('算例名称', cfg.tag);
cfg.plot     = ask_yn('运行后绘图?', sps_getp(cfg,'plot',true));
cfg.verbose  = ask_yn('打印配置摘要与诊断报告?', sps_getp(cfg,'verbose',true));
cfg.savefile = ask_str('保存结果到 .mat(留空 = 不保存)', sps_getp(cfg,'savefile',''));
end

% ======================================================================
%                            小工具
% ======================================================================
function i = ask_choice(title, names, dflt)
fprintf('\n--- %s ---\n', title);
for k = 1:numel(names)
    if k == dflt, mark = '*'; else, mark = ' '; end
    fprintf(' %s%2d) %s\n', mark, k, names{k});
end
s = input(sprintf('选择编号 [默认 %d]: ', dflt), 's');
if isempty(strtrim(s))
    i = dflt;
    return;
end
v = str2double(s);
if isnan(v) || v < 1 || v > numel(names)
    fprintf('  输入无效,使用默认 %d\n', dflt);
    i = dflt;
else
    i = round(v);
end
end

% ----------------------------------------------------------------------
function v = ask_num(prompt, dflt)
if isempty(dflt) || ~isnumeric(dflt) || ~isscalar(dflt) || ~isfinite(dflt)
    dflt = 0;
end
s = input(sprintf('%s [默认 %.6g]: ', prompt, dflt), 's');
if isempty(strtrim(s))
    v = dflt;
    return;
end
try
    vv = str2num(s); %#ok<ST2NM>
catch
    vv = [];
end
if isempty(vv) || ~isscalar(vv) || ~isfinite(vv)
    fprintf('  无法解析,保留默认 %.6g\n', dflt);
    v = dflt;
else
    v = vv;
end
end

% ----------------------------------------------------------------------
function v = ask_vec3(prompt, dflt)
if isempty(dflt)
    dflt = [0;0;0];
end
dflt = dflt(:);
if numel(dflt) ~= 3
    dflt = [0;0;0];
end
s = input(sprintf('%s [默认 %s]: ', prompt, mat2str(dflt.',4)), 's');
if isempty(strtrim(s))
    v = dflt;
    return;
end
vv = str2num(s); %#ok<ST2NM>
if numel(vv) == 3
    v = vv(:);
else
    fprintf('  需要 3 个分量,保留默认 %s\n', mat2str(dflt.',4));
    v = dflt;
end
end

% ----------------------------------------------------------------------
function s = ask_str(prompt, dflt)
if ~ischar(dflt), dflt = ''; end
t = input(sprintf('%s [默认 %s]: ', prompt, dflt), 's');
if isempty(strtrim(t))
    s = dflt;
else
    s = strtrim(t);
end
end

% ----------------------------------------------------------------------
function tf = ask_yn(prompt, dflt)
if dflt, d = 'Y/n'; else, d = 'y/N'; end
t = input(sprintf('%s [%s]: ', prompt, d), 's');
t = strtrim(lower(t));
if isempty(t)
    tf = dflt;
elseif any(strcmp(t, {'y','yes','1','是','true'}))
    tf = true;
else
    tf = false;
end
end

% ----------------------------------------------------------------------
function s = rmfield_safe(s, f)
if isfield(s, f)
    s = rmfield(s, f);
end
end
