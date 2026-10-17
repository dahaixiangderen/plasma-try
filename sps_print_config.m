function sps_print_config(cfg)
%SPS_PRINT_CONFIG  打印配置摘要(粒子、场模型、初值、回旋尺度、时间步)

c = cfg.const;
fprintf('\n==================== 配置摘要: %s ====================\n', cfg.tag);

% ---------- 粒子 ----------
KE_eV = 0.5*cfg.m*cfg.scales.v0^2/c.e;
fprintf('粒子      : %s\n', cfg.species_label);
fprintf('            q = %+.4e C,  m = %.4e kg,  相对论 = %d\n', ...
        cfg.q, cfg.m, cfg.relativistic);
fprintf('初速度    : |v0| = %.4e m/s = %.4e km/s, 动能 = %.4e eV (非相对论)\n', ...
        cfg.scales.v0, cfg.scales.v0/1e3, KE_eV);

% ---------- 位置 ----------
x0 = cfg.init.x0;
fprintf('初始位置  : [%.4e, %.4e, %.4e] m\n', x0(1), x0(2), x0(3));
fprintf('             = %s\n', pos_string(x0, c));

% ---------- 磁场 ----------
fprintf('磁场模型  : %s,  |B(x0)| = %.4e T = %.3f nT\n', ...
        cfg.Bmodel, cfg.scales.B_now, cfg.scales.B_now*1e9);
fprintf('            B0 = [%.4g, %.4g, %.4g] nT\n', cfg.B0(1)*1e9, cfg.B0(2)*1e9, cfg.B0(3)*1e9);
print_field_params('Bparams', cfg.Bparams);
if isfield(cfg.Bmod,'amp') && cfg.Bmod.amp ~= 0
    fprintf('            |B| 时间调制: amp=%.3g, omega=%.3g rad/s, phase=%.3g rad\n', ...
            cfg.Bmod.amp, cfg.Bmod.omega, cfg.Bmod.phase);
end

% ---------- 电场 ----------
[Enow, ~] = sps_fields(0, x0, cfg);
fprintf('电场模型  : %s,  |E(x0)| = %.4e V/m = %.3f mV/m\n', ...
        cfg.Emodel, norm(Enow), norm(Enow)*1e3);
fprintf('            E0 = [%.4g, %.4g, %.4g] mV/m\n', cfg.E0(1)*1e3, cfg.E0(2)*1e3, cfg.E0(3)*1e3);
print_field_params('Eparams', cfg.Eparams);

% ---------- 回旋尺度 ----------
fprintf('回旋尺度  : Omega_c = %.5g rad/s,  T_c = %.5g s,  f_c = %.5g Hz\n', ...
        cfg.scales.omega_c, cfg.scales.Tc, cfg.scales.fc);
fprintf('            r_g0 = %.5g m = %.4g km\n', cfg.scales.rg0, cfg.scales.rg0/1e3);
if isfinite(cfg.scales.Twave)
    fprintf('最快场时标: T_wave = %.5g s\n', cfg.scales.Twave);
end

% ---------- 数值 ----------
nu = cfg.numerics;
fprintf('数值      : dt = %.5g s (%.1f 步/回旋周期), 步数 = %d, 总时长 = %.5g s (%.1f T_c)\n', ...
        nu.dt, cfg.scales.Tc/nu.dt, nu.nsteps, nu.tmax, nu.tmax/cfg.scales.Tc);
fprintf('            子步 = %d, 存储间隔 = %d, 差分步长 = %.4g m\n', ...
        nu.substeps, nu.store_every, nu.diff_h);
fprintf('=========================================================\n\n');

end

% ----------------------------------------------------------------------
function print_field_params(name, p)
%PRINT_FIELD_PARAMS  打印场模型参数(跳过内部量;磁场类标量按 nT 显示)

if ~isstruct(p) || isempty(fieldnames(p))
    return;
end
skip   = {'origin','modes','phase0','phase','k'};

bnames = {'Bmag','Bx','Bn','Bg','B1','B_r0','dB_rms','moment','eq_field'};

fn = fieldnames(p);
shown = 0;
for i = 1:numel(fn)
    key = fn{i};
    if any(strcmpi(key, skip))
        continue;
    end
    val = p.(key);
    if ischar(val)
        fprintf('            %s.%s = ''%s''\n', name, key, val);
    elseif isnumeric(val) && isscalar(val) && isfinite(val)
        if any(strcmpi(key, bnames))
            fprintf('            %s.%s = %.6g nT\n', name, key, val*1e9);
        else
            fprintf('            %s.%s = %.6g\n', name, key, val);
        end
    elseif isnumeric(val) && ~isempty(val) && numel(val) <= 4
        fprintf('            %s.%s = [%s]\n', name, key, num2str(val(:).', '%g '));
    else
        continue;
    end
    shown = shown + 1;
    if shown >= 10
        fprintf('            ...\n');
        break;
    end
end

end

% ----------------------------------------------------------------------
function s = pos_string(x, c)
%POS_STRING  用 km / 地球半径 / AU 表示位置,便于阅读

r = norm(x);
if r >= 0.05*c.AU
    s = sprintf('%.4g AU', r/c.AU);
elseif r >= 0.5*c.Re
    s = sprintf('%.4g R_E (%.4g km)', r/c.Re, r/1e3);
else
    s = sprintf('%.4g km', r/1e3);
end

end
