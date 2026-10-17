function cfg = sps_prepare(cfg)
%SPS_PREPARE  校验/补全配置,预计算湍流场,并自动确定时间步长与步数
%
%   调用后 cfg 会新增/解析出以下字段:
%     cfg.q, cfg.m, cfg.species_label   粒子电荷/质量/名称
%     cfg.init.x0, cfg.init.v0          解析后的初始位置与速度
%     cfg.scales.Bref/Bmax/omega_c/Tc/fc/rg0/Twave
%     cfg.numerics.dt/tmax/nsteps/nstore/diff_h/x_max/B_max

if nargin < 1 || isempty(cfg)
    cfg = sps_default_config();
end
if ~isstruct(cfg)
    error('sps_prepare:badcfg', 'cfg 必须是结构体(可用 sps_default_config 生成)。');
end

c = sps_const();
cfg.const = c;

% -------------------- 粒子 --------------------
[q, m, label] = sps_species(cfg);
cfg.q = q;
cfg.m = m;
cfg.species_label = label;

% -------------------- 场的补充设置 --------------------
if ~isfield(cfg,'Bparams') || ~isstruct(cfg.Bparams), cfg.Bparams = struct(); end
if ~isfield(cfg,'Eparams') || ~isstruct(cfg.Eparams), cfg.Eparams = struct(); end
if ~isfield(cfg,'Bmod')    || ~isstruct(cfg.Bmod),    cfg.Bmod = struct();    end
if ~isfield(cfg,'Emod')    || ~isstruct(cfg.Emod),    cfg.Emod = struct();    end
cfg.B0 = cfg.B0(:);
cfg.E0 = cfg.E0(:);

% 偶极场:允许"赤道表面磁场 + 行星半径"代替偶极矩
if strcmpi(cfg.Bmodel, 'dipole')
    Bp = cfg.Bparams;
    if isfield(Bp,'eq_field') && ~isempty(Bp.eq_field)
        R   = sps_getp(Bp, 'R_planet', c.Re);
        dir = sps_unit(sps_getp(Bp, 'moment_dir', [0;0;-1]));
        cfg.Bparams.moment = 4*pi*R^3*Bp.eq_field/c.mu0 * dir;
    end
    if ~isfield(cfg.Bparams,'moment') || isempty(cfg.Bparams.moment)
        cfg.Bparams.moment = [0; 0; -c.m_earth];
    end
end

% 湍流场:预生成各模态(相位固定,结果可复现)
if strcmpi(cfg.Bmodel, 'turbulent')
    if ~isfield(cfg.Bparams,'modes') || isempty(cfg.Bparams.modes)
        [md, VA] = sps_make_turbulence(cfg);
        cfg.Bparams.modes = md;
        cfg.Bparams.VA    = VA;      % 供 Emodel='turbulence_inductive' 使用
    end
end

% 检查波模型与电场模型的搭配
if any(strcmpi(cfg.Emodel, {'inductive_wave','turbulence_inductive'})) && ...
        ~any(strcmpi(cfg.Bmodel, {'wave','turbulent'}))
    warning('sps_prepare:Efield', ...
        ['Emodel=''%s'' 需要 Bmodel 为 ''wave'' 或 ''turbulent'' 才有磁场扰动,', ...
         '当前 Bmodel=''%s'',感应电场将为 0。'], cfg.Emodel, cfg.Bmodel);
end

% -------------------- 初始条件(解析成直角坐标) --------------------
if ~isfield(cfg,'init') || ~isstruct(cfg.init), cfg.init = struct(); end
if ~isfield(cfg.init,'mode') || isempty(cfg.init.mode), cfg.init.mode = 'gc'; end
cfg.init.x0 = cfg.init.x0(:);
if ~isfield(cfg.init,'vdrift') || isempty(cfg.init.vdrift), cfg.init.vdrift = [0;0;0]; end
cfg.init.vdrift = cfg.init.vdrift(:);

[x0, v0] = sps_initial_state(cfg);
cfg.init.x0 = x0;
cfg.init.v0 = v0;
v0mag = norm(v0);

% -------------------- 回旋尺度 --------------------
Bmax_s = estimate_Bmax(cfg, x0);      % 采样区域内的最大 |B|(仅作参考/报告)
[~, Bx] = sps_fields(0, x0, cfg);
Bnow = norm(Bx);
if ~isfinite(Bnow)
    error('sps_prepare:sing', ...
        ['初始点处的磁场不是有限值(通常是初始位置落在偶极场/电流片的奇点上)。', ...
         '请修改 cfg.init.x0。']);
end
if Bnow > 0
    Bref = Bnow;                      % 时间步长按初始点附近的局地场确定
else
    Bref = max(Bmax_s, realmin);      % 初始点恰为零场时的退化处理
    warning('sps_prepare:zerofield', ...
        ['初始点处 |B| = 0,回旋运动退化(回旋半径无穷大),', ...
         '时间步按采样区域最大 |B| 估计,结果仅供定性参考。']);
end

omega_c = abs(q)*Bref/m;          % 回旋角频率 [rad/s]
Tc      = 2*pi/omega_c;           % 回旋周期 [s]
rg0     = v0mag/omega_c;          % 初始回旋半径 [m]

% 场中含时的最快频率(波、调制、湍流模)
oms = [ abs(sps_getp(cfg.Bparams,'omega',0));
        abs(sps_getp(cfg.Eparams,'omega',0));
        abs(sps_getp(cfg.Eparams,'Uomega',0));
        abs(sps_getp(cfg.Bmod,'omega',0));
        abs(sps_getp(cfg.Emod,'omega',0)) ];
if strcmpi(cfg.Bmodel,'turbulent')
    md = sps_getp(cfg.Bparams,'modes',[]);
    for i = 1:numel(md)
        oms(end+1,1) = abs(md(i).omega); %#ok<AGROW>
    end
end
om_max = max(oms);
if om_max > 0
    Twave = 2*pi/om_max;
else
    Twave = Inf;
end

cfg.scales = struct('Bmax',Bmax_s, 'B_now',Bnow, 'Bref',Bref, ...
                    'omega_c',omega_c, 'Tc',Tc, 'fc',1/Tc, ...
                    'rg0',rg0, 'v0',v0mag, 'Twave',Twave);

% -------------------- 时间步长与步数 --------------------
nu = cfg.numerics;
if ~isfield(nu,'nper')      || isempty(nu.nper),      nu.nper = 120;   end
if ~isfield(nu,'ngyro')     || isempty(nu.ngyro),     nu.ngyro = 30;   end
if ~isfield(nu,'substeps')  || isempty(nu.substeps),  nu.substeps = 1; end
if ~isfield(nu,'store_every')||isempty(nu.store_every),nu.store_every = 1; end
if ~isfield(nu,'adaptive')  || isempty(nu.adaptive),  nu.adaptive = true; end
if ~isfield(nu,'substeps_trigger') || isempty(nu.substeps_trigger), nu.substeps_trigger = 2; end
if ~isfield(nu,'substeps_max')||isempty(nu.substeps_max), nu.substeps_max = 200; end

if ~isempty(nu.dt) && nu.dt > 0
    dt = nu.dt;
else
    dt = min(Tc, Twave)/max(nu.nper, 10);
end

if dt > Tc/20
    warning('sps_prepare:bigdt', ...
        ['时间步 dt=%.3g s 仅能分辨 %.1f 步/回旋周期,回旋运动会有明显误差。', ...
         '建议增大 cfg.numerics.nper 或减小 dt。'], dt, Tc/dt);
end
if isfinite(Twave) && dt > Twave/20
    warning('sps_prepare:bigdtwave', ...
        'dt=%.3g s 仅能分辨 %.1f 步/波周期,波场采样不足。', dt, Twave/dt);
end

if ~isempty(nu.tmax) && nu.tmax > 0
    tmax = nu.tmax;
else
    tmax = max(nu.ngyro, 1) * Tc;
end

nsteps = max(ceil(tmax/dt), 10);
if nsteps > 5e6
    warning('sps_prepare:many', ...
        ['共 %d 步,内存/时间开销较大。可增大 cfg.numerics.store_every ', ...
         '或减小 cfg.numerics.ngyro。'], nsteps);
end

% 空间差分步长:取回旋半径的 1/100(同时受浮点精度下限约束)
if ~isempty(nu.diff_h) && nu.diff_h > 0
    diff_h = nu.diff_h;
else
    diff_h = max(rg0/100, 1);
end
diff_h = max(diff_h, max(norm(x0),1)*1e-10);

% 逃逸判据:用户显式给了就用用户的;否则按场位形与自由飞行距离自动取
if ~isempty(nu.x_max) && nu.x_max > 0
    x_max = nu.x_max;
else
    sx = max(norm(x0), 1e6);
    switch lower(cfg.Bmodel)
        case 'dipole',  sx = max(sx,  8*c.Re);
        case 'parker',  sx = max(sx,  3*c.AU);
        case 'mirror',  sx = max(sx, 50*sps_getp(cfg.Bparams,'L',1e6));
        case 'harris',  sx = max(sx, 50*sps_getp(cfg.Bparams,'L',1e6));
    end
    x_max = 50*sx;
    x_max = max(x_max, 3*v0mag*tmax + norm(x0));   % 至少允许自由飞行 3 倍行程
end

nu.dt          = dt;
nu.tmax        = tmax;
nu.nsteps      = nsteps;
nu.nstore      = floor(nsteps/nu.store_every) + 1;
nu.diff_h      = diff_h;
nu.x_max       = x_max;
nu.B_max       = 100*Bref;
nu.dt_over_Tc  = dt/Tc;
cfg.numerics   = nu;

if cfg.verbose
    sps_print_config(cfg);
end

end

% ======================================================================
function Bmax = estimate_Bmax(cfg, x0)
%ESTIMATE_BMAX  在初始点周围采样,估计模拟区域内的最大 |B|
%   用于按最强场确定时间步长(避免粒子走进强场区时步长不足)。

c = sps_const();
switch lower(cfg.Bmodel)
    case 'dipole',          scales = [0.5, 1, 2, 4, 8]*c.Re;
    case {'harris','mirror','linear_gradient'}, scales = [1e5, 1e6, 1e7];
    case 'parker',          scales = [0.5, 1, 2]*c.AU;
    case 'turbulent',       scales = [1e5, 1e6, 1e7];
    otherwise,              scales = [1e6, 1e7, 1e8];
end

pts = x0(:);
for s = scales
    for k = 1:3
        e = zeros(3,1);
        e(k) = s;
        pts = [pts, x0(:)+e, x0(:)-e]; %#ok<AGROW>
    end
end

Bmax = norm(cfg.B0(:));
for i = 1:size(pts,2)
    if strcmpi(cfg.Bmodel,'dipole')
        o = sps_getp(cfg.Bparams,'origin',[0;0;0]);
        if norm(pts(:,i) - o(:)) < 0.2*c.Re
            continue;                       % 太靠近奇点,跳过
        end
    end
    [~, Bi] = sps_fields(0, pts(:,i), cfg);
    bi = norm(Bi);
    if isfinite(bi) && bi > 0 && bi < 1e2
        Bmax = max(Bmax, bi);
    end
end

if ~isfinite(Bmax) || Bmax <= 0
    [~, Bi] = sps_fields(0, x0, cfg);
    Bmax = max(norm(Bi), norm(cfg.B0(:)));
end
if ~isfinite(Bmax) || Bmax <= 0
    Bmax = 1e-6;                            % 兜底
end

end
