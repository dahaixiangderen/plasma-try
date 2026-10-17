function cfg = sps_default_config()
%SPS_DEFAULT_CONFIG  返回一份完整的默认配置结构体
%
%   所有可选项都在这里列出;任何预设/菜单/手动脚本只需覆盖需要的字段。
%   完整选项说明见 README.md 与 sps_fields.m 的注释。
%
%   单位约定(全部 SI):
%     位置 m, 速度 m/s, 时间 s, 磁场 T, 电场 V/m, 投掷角/相位 deg

cfg = struct();

cfg.tag          = 'default';     % 算例名(用于保存文件名/图标题)
cfg.species      = 'proton';      % 'proton'|'electron'|'alpha'|'O+'|'custom'
cfg.q            = [];            % species='custom' 时的电荷 [C](带符号)
cfg.m            = [];            % species='custom' 时的质量 [kg]
cfg.relativistic = false;         % true: 使用相对论 Boris 推进(高能电子等)

% ======================= 磁场模型 =======================
% Bmodel 可选:
%   'uniform'         均匀恒定磁场 B0
%   'linear_gradient' |B| 沿给定方向线性变化(梯度-B 漂移)
%   'dipole'          地球/行星偶极场(弯曲场线:梯度-B + 曲率漂移)
%   'harris'          Harris 电流片(磁尾重联位形,Bx=tanh(z/L))
%   'mirror'          轴对称磁镜(无散度,演示磁镜捕获/损失锥)
%   'wave'            背景场 + 圆/线偏振平面波(哨声波、EMIC 等)
%   'parker'          Parker 螺旋太阳风磁场
%   'turbulent'       背景场 + 多模 Alfven 湍流(随机相位,可复现)
%   'custom'          用户自定义函数 cfg.Bfun(t,x)
cfg.Bmodel  = 'uniform';
cfg.B0      = [0; 0; 50e-9];      % 背景磁场 [T], 3x1
cfg.Bparams = struct();           % 模型专属参数(见 sps_fields.m)
cfg.Bmod    = struct('amp',0, 'omega',0, 'phase',0);  % |B| 整体时间调制
%   B(t) = B_model(t,x) * (1 + amp*sin(omega*t + phase))
%   注意: 这是"运动学"调制,不含 Faraday 定律要求的感应电场,
%   因此它单独使用时不是自洽的电磁场(磁矩会随 1/|B| 变化)。
%   需要自洽的时变场请用 Bmodel='wave' + Emodel='inductive_wave'。

% ======================= 电场模型 =======================
% Emodel 可选:
%   'zero'                   E = 0
%   'uniform'                均匀静电场 E0
%   'time_harmonic'          E0 + E1*cos(omega*t+phase) 沿 Edir
%   'circular'               垂直于 B 的圆偏振电场(回旋共振加热)
%   'convective'             对流电场 E = -U x B (太阳风/磁层对流)
%   'inductive_wave'         由 Faraday 定律给出的波感应电场
%   'parallel_wave'          沿 B 的平行电场(朗道加速/平行电场)
%   'turbulence_inductive'   湍流模对应的感应电场(与 Bmodel='turbulent' 配套)
%   'custom'                 用户自定义函数 cfg.Efun(t,x)
cfg.Emodel  = 'zero';
cfg.E0      = [0; 0; 0];          % 静电场 [V/m], 3x1
cfg.Eparams = struct();           % 模型专属参数(见 sps_fields.m)
cfg.Emod    = struct('amp',0, 'omega',0, 'phase',0);  % E 整体时间调制

% ======================= 初始条件 =======================
% init.mode = 'gc'  : 用 速率/投掷角/回旋相位 给初速度(推荐,物理直观)
%           = 'cart': 直接给直角坐标位置与速度
cfg.init = struct( ...
    'mode',      'gc',   ...   % 'gc' | 'cart'
    'x0',        [0;0;0], ...  % 初始位置 [m]
    'v0',        [0;0;0], ...  % mode='cart' 时的初始速度 [m/s]
    'speed',     1e5,   ...    % 速率 [m/s]
    'alpha',     90,    ...    % 投掷角 alpha [deg] (v_par = v cos, v_perp = v sin)
    'gyrophase', 0,     ...    % 回旋相位 [deg]
    'vpar',      NaN,   ...    % 若为有限值则直接用 v_par,忽略 speed/alpha
    'vperp',     NaN,   ...    % 若为有限值则直接用 v_perp
    'vdrift',    [0;0;0]);     % 附加整体漂移速度 [m/s]

% ======================= 数值设置 =======================
cfg.numerics = struct( ...
    'nper',        120, ...  % 每回旋周期步数(用于自动确定 dt)
    'dt',          [],  ...  % 手动步长 [s](非空则覆盖 nper)
    'tmax',        [],  ...  % 手动总时长 [s](非空则覆盖 ngyro)
    'ngyro',       30,  ...  % 模拟多少个初始回旋周期
    'substeps',    1,   ...  % 每个 dt 内的 Boris 子步数(最小值)
    'adaptive',    true, ...  % 若局地 |B| 明显变强则自动增加子步(安全网)
    'substeps_trigger', 2, ...% 局地回旋周期低于 dt*trigger 时才细分
    'substeps_max',200, ...  % 自适应子步上限
    'diff_h',      [],  ...  % 空间差分步长 [m](空=按回旋半径自动取)
    'store_every', 1,   ...  % 每多少步存一次(节省内存)
    'x_max',       []);      % 位置上限 [m],超出即判定粒子逃逸并提前结束

cfg.seed     = 20240501;  % 随机种子(湍流场相位)
cfg.verbose  = true;      % 打印配置摘要与诊断报告
cfg.plot     = true;      % 是否绘图
cfg.savefile = '';        % 非空则把结果保存到该 .mat 文件

end
