function out = sps_analysis(out)
%SPS_ANALYSIS  对模拟结果做定量诊断:守恒量、导心漂移、回旋/反弹周期
%
%   结果写入 out.summary,并在 cfg.verbose 为真时打印诊断报告。
%
%   诊断内容:
%     1) 回旋尺度(理论)与实测回旋周期
%     2) 磁矩 mu 与动能 KE 的相对漂移(数值守恒性检验)
%     3) 各导心漂移分量(解析)与模拟得到的导心位移(实测)对比
%     4) 平行速度变号检测 -> 磁镜反弹周期/捕获判据

cfg = out.cfg;
c   = cfg.const;
s   = struct();

nt = numel(out.t);
s.npoints = nt;
s.omega_c = cfg.scales.omega_c;
s.Tc      = cfg.scales.Tc;
s.rg0     = cfg.scales.rg0;

% ---------------- 守恒性 ----------------
s.KE0     = out.KE(1);
s.KE0_eV  = out.KE(1)/c.e;
s.dKE_rel = max(abs(out.KE - out.KE(1)))/max(abs(out.KE(1)), realmin);
if all(isfinite(out.mu))
    s.mu0     = out.mu(1);
    s.dmu_rel = max(abs(out.mu - out.mu(1)))/max(abs(out.mu(1)), realmin);
else
    s.mu0 = NaN;
    s.dmu_rel = NaN;
end

% 回旋平均后的磁矩:去掉 O(r_g*|grad B|/B) 量级的高频摆动,
% 剩下的才是真正的第一绝热不变量(应在漂移时间尺度上守恒)。
% 窗口取 2 个回旋周期,并且只统计窗口完整的内部点(边界处平均不完整)。
nsm = max(5, round(2*cfg.numerics.nper/max(1, round(cfg.numerics.store_every))));
if nt >= nsm + 10 && all(isfinite(out.mu))
    mus  = movmean(out.mu, nsm);
    i0   = floor(nsm/2);
    core = mus(i0+1 : nt-i0);
    s.mu_smooth       = mus;
    s.mu_core         = core;
    s.mu_smooth_mean  = mean(core);
    s.dmu_smooth_rel  = max(abs(core - core(1)))/max(abs(core(1)), realmin);
else
    s.mu_smooth      = out.mu;
    s.mu_core        = out.mu;
    s.mu_smooth_mean = mean(out.mu(isfinite(out.mu)));
    s.dmu_smooth_rel = s.dmu_rel;
end

% ---------------- 实测回旋周期(垂直速度在固定方向上的投影过零) ----------------
% 只有当磁场在整条轨道上"方向基本不变、强度基本不变"时,固定方向的投影过零
% 才对应以初始点 Omega_c 定义的回旋周期。否则(磁镜里 |B| 变化 2 倍、
% 强弯曲场线里方向大幅转动)这个测量没有可比性,直接跳过。
s.Tc_measured = NaN;
s.Tc_measured_valid = false;
if nt >= 5
    [~, B1] = sps_fields(out.t(1), out.x(:,1), cfg);
    bh0 = sps_unit(B1);
    bn_all = sqrt(sum(out.B.^2, 1));
    if norm(bh0) > 0 && all(bn_all > 0)
        rot    = acosd(min(max((out.B ./ bn_all).' * bh0, -1), 1));
        Bvar   = max(bn_all)/min(bn_all);
    else
        rot  = 180;
        Bvar = Inf;
    end
    if norm(bh0) > 0 && max(rot) < 20 && Bvar < 1.1
        [e1, ~] = sps_perp_basis(bh0, bh0);
        proj = e1.' * out.v;
        sg = sign(proj);
        sg(sg == 0) = 1;
        nc = sum(abs(diff(sg)) > 0);
        if nc >= 4
            s.Tc_measured = 2*(out.t(end) - out.t(1))/nc;
            s.Tc_measured_valid = true;
        end
    end
end

% ---------------- 导心漂移:解析 vs 实测 ----------------
if nt >= 2
    dta = out.t(end) - out.t(1);
    s.vgc_measured = (out.gc(:,end) - out.gc(:,1))/dta;   % 导心位移/时间
    s.vx_measured  = (out.x(:,end)  - out.x(:,1)) /dta;   % 真实位移/时间
    d1 = sps_drifts(out.t(1),    out.x(:,1),   out.v(:,1),   cfg);
    d2 = sps_drifts(out.t(end),  out.x(:,end), out.v(:,end), cfg);
    s.drift_initial = d1;
    s.drift_final   = d2;
    s.vdrift_analytic = 0.5*(d1.vtotal + d2.vtotal);
    s.drift_speed_measured = norm(s.vgc_measured);
    s.drift_speed_analytic = norm(s.vdrift_analytic);

    % 绕 z 轴的漂移周期(仅在导心确有方位角运动时有意义)
    rperp = norm([out.gc(1,1); out.gc(2,1)]);
    if rperp > 0
        ephi = [-out.gc(2,1); out.gc(1,1); 0]/rperp;
        vphi = dot(s.vgc_measured, ephi);
        if abs(vphi) > 1e-9
            s.drift_period_azimuthal = 2*pi*rperp/abs(vphi);
            s.drift_direction = sign(vphi);
        end
    end
else
    s.vgc_measured = [NaN;NaN;NaN];
    s.vx_measured  = [NaN;NaN;NaN];
    s.drift_initial = struct();
    s.drift_final   = struct();
    s.vdrift_analytic = [NaN;NaN;NaN];
    s.drift_speed_measured = NaN;
    s.drift_speed_analytic = NaN;
end

% ---------------- 磁镜反弹 / 损失锥 ----------------
s.bounce_period = NaN;
s.n_bounce = 0;
vpars = cosd(out.pitch) .* out.gamma;              % 与 v_par 同号
vscale = max(abs(vpars));
if vscale > 0
    sgn = zeros(size(vpars));
    sgn(vpars >  0.01*vscale) =  1;
    sgn(vpars < -0.01*vscale) = -1;
    sel = find(sgn ~= 0);
    if numel(sel) > 2
        sz = sgn(sel);
        tz = out.t(sel);
        ci = find(diff(sz) ~= 0);
        if numel(ci) >= 2
            tc = tz(ci);
            s.bounce_period = 2*median(diff(tc));
            s.n_bounce = numel(ci) - 1;
        end
    end
end
s.pitch_min = min(out.pitch);
s.pitch_max = max(out.pitch);
s.pitch_spread = s.pitch_max - s.pitch_min;

% ---------------- 能量/速度极值 ----------------
s.vmax = max(sqrt(sum(out.v.^2, 1)));
s.vmin = min(sqrt(sum(out.v.^2, 1)));
s.Bmag_min = min(out.Bmag);
s.Bmag_max = max(out.Bmag);
s.Emag_max = max(out.Emag);
s.x_span = max(out.x, [], 2) - min(out.x, [], 2);
s.tmax_sim = out.t(end) - out.t(1);

out.summary = s;

% ---------------- 报告 ----------------
if cfg.verbose
    fprintf('\n--------------------- 诊断报告: %s ---------------------\n', cfg.tag);
    fprintf('模拟: %d 个存储点, 时长 %.4g s (%.2f 个回旋周期), 用时 %.2f s\n', ...
            nt, s.tmax_sim, s.tmax_sim/s.Tc, out.walltime);
    if out.escaped
        fprintf('>> 轨道提前终止: %s\n', out.esc_reason);
    end
    fprintf('\n[回旋运动]\n');
    fprintf('  理论 T_c = %.5g s (Omega_c = %.5g rad/s),  r_g0 = %.4g km\n', ...
            s.Tc, s.omega_c, s.rg0/1e3);
    if isfinite(s.Tc_measured)
        fprintf('  实测 T_c = %.5g s (相对误差 %.3g)\n', ...
                s.Tc_measured, abs(s.Tc_measured - s.Tc)/s.Tc);
    elseif ~s.Tc_measured_valid
        fprintf('  实测 T_c : 轨道上磁场方向或强度变化较大,该项不适用(未测量)\n');
    end
    fprintf('  模拟区域 |B| 范围: %.4g ~ %.4g nT ; |E| 最大 %.4g mV/m\n', ...
            s.Bmag_min*1e9, s.Bmag_max*1e9, s.Emag_max*1e3);

    fprintf('\n[守恒量]\n');
    fprintf('  动能 KE0 = %.4g J = %.4g eV,  |dKE|/KE = %.3e\n', ...
            s.KE0, s.KE0_eV, s.dKE_rel);
    fprintf('  磁矩 mu0 = %.4g J/T,  |dmu|/mu = %.3e (瞬时值,含回旋摆动)\n', ...
            s.mu0, s.dmu_rel);
    fprintf('  回旋平均后        ,  |dmu|/mu = %.3e (真正的绝热不变量,均值 %.5g J/T)\n', ...
            s.dmu_smooth_rel, s.mu_smooth_mean);
    if s.dKE_rel > 1e-3 && s.Emag_max > 0
        fprintf('  (注: 存在电场做功,动能本就不守恒;应重点看 mu 的守恒性)\n');
    end

    fprintf('\n[导心漂移](单位 m/s)\n');
    if isfield(s.drift_initial,'vExB')
        d = s.drift_initial;
        fprintf('  解析(初始点): E x B = %.4g | grad-B = %.4g | 曲率 = %.4g | 平行 = %.4g | 极化 = %.3g\n', ...
                norm(d.vExB), norm(d.vgradB), norm(d.vcurv), norm(d.vparB), norm(d.vpol));
        fprintf('                 合成 v_d = %.4g\n', norm(d.vtotal));
    end
    fprintf('  实测(导心位移/时间): |v_gc| = %.4g,  |v_particle| = %.4g\n', ...
            s.drift_speed_measured, norm(s.vx_measured));
    if isfinite(s.drift_speed_measured) && isfinite(s.drift_speed_analytic) && ...
            s.drift_speed_analytic > 0
        rr = s.drift_speed_measured/s.drift_speed_analytic;
        if rr < 0.33 || rr > 3
            fprintf(['  (注: 实测/解析 = %.2g,偏差较大。通常说明回旋半径与场的变化尺度\n', ...
                     '   可比(如薄电流片中心/零场点附近),导心展开的精度有限。)\n'], rr);
        end
    end
    if isfield(s,'drift_period_azimuthal')
        fprintf('  绕轴漂移周期 = %.4g s (%.3g 分钟), 方向 = %+d\n', ...
                s.drift_period_azimuthal, s.drift_period_azimuthal/60, s.drift_direction);
    end

    fprintf('\n[投掷角/磁镜]\n');
    fprintf('  投掷角范围: %.3g ~ %.3g deg (变化 %.3g deg)\n', ...
            s.pitch_min, s.pitch_max, s.pitch_spread);
    if isfinite(s.bounce_period)
        fprintf('  检测到 %d 次平行速度变号,反弹周期约 %.4g s\n', s.n_bounce, s.bounce_period);
        if strcmpi(cfg.Bmodel,'mirror')
            fprintf('  -> 粒子被磁镜捕获(反弹运动)\n');
        end
        fprintf(['  (注: 捕获/反弹轨道的端点位移会因往返而相消,\n', ...
                 '   上面"实测导心漂移"对这类轨道不具意义)\n']);
    elseif strcmpi(cfg.Bmodel,'mirror')
        fprintf('  -> 未检测到反弹:投掷角可能落在损失锥内(粒子沿场线逃逸)\n');
    end
    fprintf('-----------------------------------------------------------\n\n');
end

end
