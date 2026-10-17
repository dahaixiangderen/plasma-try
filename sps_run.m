function out = sps_run(cfg, doPlot)
%SPS_RUN  主求解器:用 Boris 推进器积分单粒子轨道并做诊断分析
%
%   out = SPS_RUN(cfg)
%   out = SPS_RUN(cfg, doPlot)   % doPlot 覆盖 cfg.plot
%
%   输出 out 的字段:
%     out.cfg            最终使用的配置
%     out.t              时间 [s]        (1 x N)
%     out.x, out.v       位置 [m]、速度 [m/s] (3 x N)
%     out.E, out.B       轨道上的电场 [V/m]、磁场 [T] (3 x N)
%     out.Bmag, out.Emag |B| [T]、|E| [V/m]
%     out.KE             动能 [J]; out.mu 第一绝热不变量 [J/T]
%     out.pitch          投掷角 [deg]; out.gamma 相对论因子
%     out.gc             瞬时导心位置 [m] (3 x N)
%     out.escaped        是否提前终止; out.esc_reason 原因
%     out.summary        诊断汇总结构体(sps_analysis 生成)
%
%   例:
%     cfg = sps_presets('radiation_belt'); out = sps_run(cfg);

if nargin < 1 || isempty(cfg)
    cfg = sps_default_config();
end
if nargin < 2 || isempty(doPlot)
    doPlot = sps_getp(cfg, 'plot', true);
end

cfg = sps_prepare(cfg);

q   = cfg.q;
m   = cfg.m;
dt  = cfg.numerics.dt;
ns  = cfg.numerics.nsteps;
se  = max(1, round(cfg.numerics.store_every));
rel = cfg.relativistic;
sub = cfg.numerics.substeps;
xmax = cfg.numerics.x_max;
Bmax = cfg.numerics.B_max;
cc  = cfg.const.c_light;

[x, v] = sps_initial_state(cfg);
t = 0;

N = ceil(ns/se) + 3;                        % 存储点数上界(含末尾补点)
T  = zeros(1, N);
X  = zeros(3, N);
V  = zeros(3, N);
Bs = zeros(3, N);
Es = zeros(3, N);
GC = zeros(3, N);
Bm = zeros(1, N);
Em = zeros(1, N);
KE = zeros(1, N);
MU = zeros(1, N);
PA = zeros(1, N);
GA = ones(1, N);

k = 1;
store = true;
escaped = false;
esc_reason = '';
nsub_used = 0;
nsub_max  = 0;
t0 = tic;

for n = 0:ns-1

    th = t + dt/2;                          % 时间居中半步
    [Eh, Bh] = sps_fields(th, x, cfg);
    Bhn = norm(Bh);

    if store && k <= N
        T(k)  = th;
        X(:,k) = x;
        V(:,k) = v;
        Es(:,k) = Eh;
        Bs(:,k) = Bh;
        Bm(k) = Bhn;
        Em(k) = norm(Eh);

        if rel
            g = 1/sqrt(max(1 - sum(v.^2)/cc^2, 1e-12));
        else
            g = 1;
        end
        GA(k) = g;
        pmag = g*m*norm(v);
        vpar = dot(v, sps_unit(Bh));
        pperp2 = max(pmag^2 - (g*m*vpar)^2, 0);
        if rel
            KE(k) = (g - 1)*m*cc^2;          % 相对论动能
        else
            KE(k) = 0.5*m*sum(v.^2);         % 非相对论动能
        end
        if Bhn > 1e-14
            MU(k) = pperp2/(2*m*Bhn);
            PA(k) = atan2d(sqrt(max(sum(v.^2) - vpar^2, 0)), vpar);
            GC(:,k) = x + (m/(q*Bhn^2))*cross(v, Bh);
        else
            MU(k) = NaN;
            PA(k) = NaN;
            GC(:,k) = x;
        end
        k = k + 1;
        store = false;
    end

    if Bhn > Bmax
        escaped = true;
        esc_reason = sprintf('|B| = %.3g T 超过阈值 %.3g T(粒子进入强场区)', Bhn, Bmax);
        break;
    end

    % 自适应子步:局地 |B| 明显强于参考场时才细分(安全网)。
    % 阈值不能太小:子步数随粒子位置变化会让积分格式在轨道上不对称,
    % 从而引入虚假的磁矩漂移(数值实验:阈值取 1 时偶极场 |dmu|/mu 从
    % 2e-3 恶化到 7e-2)。默认阈值 2,即局地场超过参考场 2 倍才细分。
    nsub = sub;
    if cfg.numerics.adaptive && Bhn > 0
        Tloc = 2*pi*m/(abs(q)*Bhn);
        need = dt/(Tloc/max(cfg.numerics.nper,1));
        if need > cfg.numerics.substeps_trigger
            nsub = min(ceil(need), cfg.numerics.substeps_max);
        end
    end
    nsub_used = nsub_used + nsub;
    nsub_max  = max(nsub_max, nsub);

    [x, v] = sps_boris(x, v, q, m, Eh, Bh, dt, rel, nsub);
    t = t + dt;

    if ~all(isfinite(x)) || ~all(isfinite(v))
        escaped = true;
        esc_reason = '出现非有限值(数值失稳,通常是步长不足或场奇点)';
        break;
    end
    if norm(x) > xmax
        escaped = true;
        esc_reason = sprintf('粒子逃逸:|x| = %.4g m 超过阈值 %.4g m', norm(x), xmax);
        break;
    end

    if mod(n+1, se) == 0
        store = true;
    end

    if cfg.verbose && mod(n+1, max(1, floor(ns/10))) == 0
        fprintf('  进度 %3.0f%%   t = %.4g s\n', 100*(n+1)/ns, t);
    end
end

% ---- 末尾补一个存储点(未逃逸时) ----
if ~escaped && k <= N
    [Ee, Be] = sps_fields(t, x, cfg);
    Bhn = norm(Be);
    T(k) = t; X(:,k) = x; V(:,k) = v; Es(:,k) = Ee; Bs(:,k) = Be;
    Bm(k) = Bhn; Em(k) = norm(Ee);
    if rel
        g = 1/sqrt(max(1 - sum(v.^2)/cc^2, 1e-12));
    else
        g = 1;
    end
    GA(k) = g;
    pmag = g*m*norm(v);
    vpar = dot(v, sps_unit(Be));
    pperp2 = max(pmag^2 - (g*m*vpar)^2, 0);
    if rel
        KE(k) = (g - 1)*m*cc^2;
    else
        KE(k) = 0.5*m*sum(v.^2);
    end
    if Bhn > 1e-14
        MU(k) = pperp2/(2*m*Bhn);
        PA(k) = atan2d(sqrt(max(sum(v.^2)-vpar^2,0)), vpar);
        GC(:,k) = x + (m/(q*Bhn^2))*cross(v, Be);
    else
        MU(k) = NaN; PA(k) = NaN; GC(:,k) = x;
    end
    k = k + 1;
end

out = struct();
out.cfg   = cfg;
out.t     = T(1:k-1);
out.x     = X(:, 1:k-1);
out.v     = V(:, 1:k-1);
out.E     = Es(:, 1:k-1);
out.B     = Bs(:, 1:k-1);
out.gc    = GC(:, 1:k-1);
out.Bmag  = Bm(1:k-1);
out.Emag  = Em(1:k-1);
out.KE    = KE(1:k-1);
out.mu    = MU(1:k-1);
out.pitch = PA(1:k-1);
out.gamma = GA(1:k-1);
out.escaped    = escaped;
out.esc_reason = esc_reason;
out.walltime   = toc(t0);
out.nsub_total = nsub_used;
out.nsub_max   = nsub_max;

out = sps_analysis(out);

if doPlot
    sps_plots(out);
end

if isfield(cfg,'savefile') && ~isempty(cfg.savefile)
    save(cfg.savefile, 'out');
    if cfg.verbose
        fprintf('结果已保存到 %s\n', cfg.savefile);
    end
end

end
