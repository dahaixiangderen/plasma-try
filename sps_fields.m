function [E, B, aux] = sps_fields(t, x, cfg)
%SPS_FIELDS  计算 t 时刻、位置 x 处的电场 E [V/m] 和磁场 B [T]
%
%   [E, B]      = SPS_FIELDS(t, x, cfg)
%   [E, B, aux] = SPS_FIELDS(t, x, cfg)
%
%   输入:
%     t   : 时间 [s](标量)
%     x   : 位置 [m](3x1 或 1x3)
%     cfg : 由 sps_default_config / 预设 / 菜单生成的配置结构体
%
%   输出:
%     E   : 电场 [V/m] (3x1)
%     B   : 磁场 [T]   (3x1)
%     aux : 诊断附加量,含
%           aux.Bmag   |B|
%           aux.bhat   磁场单位方向
%           aux.Bwave  波的磁场扰动量(用于感应电场)
%           aux.xs     相对场源原点的坐标
%
%   ------------------------------------------------------------------
%   磁场模型 (cfg.Bmodel) 及参数 (cfg.Bparams.*)
%   ------------------------------------------------------------------
%   'uniform'         B = B0
%   'linear_gradient' |B| = Bmag*(1 + grad_amp*(gdir.xs)/L),方向 Bdir
%                    参数: Bdir, Bmag, grad_dir, grad_amp, L, origin
%   'dipole'          偶极场 B = mu0/4pi*(3(m.rh)rh - m)/r^3
%                    参数: moment(3x1, A*m^2), origin
%   'harris'          电流片: B = Bx*tanh(s)*that + Bg*ghat + Bn*nhat
%                    参数: Bx, Bg, Bn, L(半厚度), nhat(法向), tangent
%   'mirror'          无散度轴对称磁镜
%                    B = Bmag*(1+(zl/L)^2)*ez - Bmag*rho*zl/L^2 * e_rho
%                    参数: Bmag, L, axis, origin
%   'wave'            背景场 + 平面波
%                    参数: k(3x1), omega, B1, polarization('RH'|'LH'|'linear'),
%                          phase0, e1(可选,自定义垂直基)
%   'parker'          Parker 螺旋场(以 origin 为太阳中心)
%                    参数: B_r0(r0 处径向场), r0, Vsw, OmegaSun, axis
%   'turbulent'       背景场 + 多模 Alfven 湍流(cfg.Bparams.modes)
%   'custom'          B = cfg.Bfun(t,x)
%
%   所有模型之后再乘以整体时间调制 (1 + Bmod.amp*sin(Bmod.omega*t+phase))
%   注意: Bmod 是运动学调制,不含 Faraday 感应电场,单独使用时不构成自洽
%         电磁场(磁矩会随 1/|B| 变化)。自洽的时变场请用
%         Bmodel='wave' + Emodel='inductive_wave'。
%
%   ------------------------------------------------------------------
%   电场模型 (cfg.Emodel) 及参数 (cfg.Eparams.*)
%   ------------------------------------------------------------------
%   'zero'                E = 0
%   'uniform'             E = E0
%   'time_harmonic'       E = E0 + E1*cos(omega*t+phase)*Edir
%   'circular'            E = E1*(e1*cos(omega*t+phase) + e2*sin(...)),
%                         e1,e2 为垂直 B0 的右手基
%   'convective'          E = -U x B         参数: U(3x1), Uamp, Uomega
%   'inductive_wave'      E = -(omega/|k|) khat x B_wave  (Faraday 定律)
%   'parallel_wave'       E = E1*bhat*cos(k.xs - omega*t)  (需 B 已算出)
%   'turbulence_inductive'各湍流模的感应电场之和
%   'custom'              E = cfg.Efun(t,x)
%
%   最后再乘以 E 的整体时间调制 (1 + Emod.amp*sin(Emod.omega*t+phase))
%
%   说明: 本函数每次调用都要算一遍场,是仿真的热点。若自定义模型很重,
%         可在 sps_run 中减少场求值次数或降低步数。

x  = x(:);
mu0 = 4*pi*1e-7;
B0  = cfg.B0(:);
b0  = sps_unit(B0);
Bp  = cfg.Bparams;
Ep  = cfg.Eparams;

origin = sps_getp(Bp, 'origin', [0;0;0]);
origin = origin(:);
xs = x - origin;                 % 相对场源坐标

Bwave = zeros(3,1);              % 波的磁场扰动(供感应电场使用)

% =====================================================================
%                              磁场
% =====================================================================
switch lower(cfg.Bmodel)

    case 'uniform'
        B = B0;

    case 'linear_gradient'
        gdir = sps_unit(sps_getp(Bp, 'grad_dir', [1;0;0]));
        L    = sps_getp(Bp, 'L',    1e6);
        amp  = sps_getp(Bp, 'grad_amp', 0.5);
        Bdir = sps_unit(sps_getp(Bp, 'Bdir', b0));
        Bmg  = sps_getp(Bp, 'Bmag', norm(B0));
        % 把梯度方向投影到垂直磁场方向,保证 div B = 0(纯平行梯度无法自洽)
        gdir = gdir - dot(gdir, Bdir)*Bdir;
        if norm(gdir) < 1e-12
            warning('sps_fields:graddir', ...
                'grad_dir 与 Bdir 平行,已退化为均匀场(平行梯度无法满足 div B = 0)。');
            gdir = zeros(3,1);
        else
            gdir = gdir/norm(gdir);
        end
        s    = dot(gdir, xs) / L;
        B    = Bdir * Bmg * (1 + amp*s);

    case 'dipole'
        mvec = sps_getp(Bp, 'moment', [0; 0; -8.07e22]);
        mvec = mvec(:);
        rn   = max(norm(xs), eps);
        rh   = xs / rn;
        B    = mu0/(4*pi) * (3*dot(mvec,rh)*rh - mvec) / rn^3;

    case 'harris'
        L     = sps_getp(Bp, 'L', 1e6);
        Bx    = sps_getp(Bp, 'Bx', 20e-9);
        Bg    = sps_getp(Bp, 'Bg', 0);
        Bn    = sps_getp(Bp, 'Bn', 0);
        nhat  = sps_unit(sps_getp(Bp, 'nhat', [0;0;1]));
        that  = sps_getp(Bp, 'tangent', [1;0;0]);
        that  = that(:);
        that  = sps_unit(that - dot(that,nhat)*nhat);   % 正交化
        ghat  = cross(nhat, that);
        s     = dot(xs, nhat) / L;
        B     = Bx*tanh(s)*that + Bg*ghat + Bn*nhat;

    case 'mirror'
        % 有限长度轴对称磁镜(磁瓶),场强在 |z| = L 处达到最大 R*Bmag:
        %   B_z = Bmag*[1 + (R-1)*sin^2(pi*z/(2L))]   (|z| <= L)
        %   B_z = Bmag*R                              (|z| >  L)
        %   B_r = -(rho/2) dB_z/dz                     (保证 div B = 0)
        % 镜比 R = B_max/B_min,损失锥角 alpha_loss = asin(sqrt(1/R))
        ez   = sps_unit(sps_getp(Bp, 'axis', [0;0;1]));
        L    = sps_getp(Bp, 'L', 1e6);
        Bmg  = sps_getp(Bp, 'Bmag', norm(B0));
        R    = max(sps_getp(Bp, 'ratio', 4), 1);
        zl   = dot(xs, ez);
        rvec = xs - zl*ez;
        rho  = norm(rvec);
        if abs(zl) <= L
            f  = 1 + (R-1)*sin(pi*zl/(2*L))^2;
            fp = (R-1)*(pi/(2*L))*sin(pi*zl/L);
        else
            f  = R;
            fp = 0;
        end
        Bz   = Bmg*f;
        Brho = -Bmg*fp*rho/2;
        if rho > 0
            B = Bz*ez + Brho*(rvec/rho);
        else
            B = Bz*ez;
        end

    case 'wave'
        kv    = sps_getp(Bp, 'k', [0;0;0]);  kv = kv(:);
        om    = sps_getp(Bp, 'omega', 0);
        B1    = sps_getp(Bp, 'B1', 0);
        pol   = sps_getp(Bp, 'polarization', 'RH');
        ph0   = sps_getp(Bp, 'phase0', 0);
        e1u   = sps_getp(Bp, 'e1', []);
        if isempty(e1u)
            [e1, e2] = sps_perp_basis(kv, b0);
        else
            kh = sps_unit(kv);
            if norm(kh) == 0, kh = sps_unit(b0); end
            if norm(kh) == 0, kh = [0;0;1];      end
            e1u = e1u(:);
            e1  = sps_unit(e1u - dot(e1u,kh)*kh);   % 投影到垂直 k 的平面
            if norm(e1) == 0
                [e1, e2] = sps_perp_basis(kv, b0);
            else
                e2 = cross(kh, e1);                 % 保证 (e1,e2,kh) 右手正交
            end
        end
        ph = dot(kv, xs) - om*t + ph0;
        switch lower(pol)
            case {'rh','right','rhcp'}
                Bwave = B1*(e1*cos(ph) + e2*sin(ph));
            case {'lh','left','lhcp'}
                Bwave = B1*(e1*cos(ph) - e2*sin(ph));
            case {'linear','lin'}
                Bwave = B1*e1*cos(ph);
            otherwise
                error('sps_fields:polarization', ...
                      '未知极化 "%s" (可选 RH/LH/linear)', pol);
        end
        B = B0 + Bwave;

    case 'parker'
        ez   = sps_unit(sps_getp(Bp, 'axis', [0;0;1]));
        r0   = sps_getp(Bp, 'r0', 1.495978707e11);
        Br0  = sps_getp(Bp, 'B_r0', 3e-9);
        Vsw  = sps_getp(Bp, 'Vsw', 4e5);
        Om   = sps_getp(Bp, 'OmegaSun', 2.865e-6);
        rn   = max(norm(xs), eps);
        rh   = xs / rn;
        sth  = norm(rh - dot(rh,ez)*ez);        % sin(theta)
        Br   = Br0 * (r0/rn)^2;
        Bph  = -Br * Om * (rn - r0) * sth / max(Vsw, eps);
        phih = sps_unit(cross(ez, rh));
        B    = Br*rh + Bph*phih;

    case 'turbulent'
        modes = sps_getp(Bp, 'modes', []);
        acc   = zeros(3,1);
        for i = 1:numel(modes)
            mm = modes(i);
            ph = dot(mm.k, xs) - mm.omega*t + mm.phase;
            acc = acc + mm.dB*(mm.e1*cos(ph) + mm.e2*sin(ph));
        end
        Bwave = acc;
        B     = B0 + acc;

    case 'custom'
        if ~isfield(cfg,'Bfun') || isempty(cfg.Bfun)
            error('sps_fields:customB', 'Bmodel=''custom'' 需要 cfg.Bfun = @(t,x) ...');
        end
        B = cfg.Bfun(t, x);
        B = B(:);
        Bwave = B - B0;

    otherwise
        error('sps_fields:Bmodel', '未知磁场模型 "%s"', cfg.Bmodel);
end

% ---- |B| 的整体时间调制(磁压缩/呼吸模) ----
if isfield(cfg,'Bmod') && isstruct(cfg.Bmod)
    a = sps_getp(cfg.Bmod, 'amp', 0);
    if a ~= 0
        B = B * (1 + a*sin(sps_getp(cfg.Bmod,'omega',0)*t + sps_getp(cfg.Bmod,'phase',0)));
    end
end

Bn = norm(B);
bh = sps_unit(B);

% =====================================================================
%                              电场
% =====================================================================
switch lower(cfg.Emodel)

    case 'zero'
        E = zeros(3,1);

    case 'uniform'
        E = cfg.E0(:);

    case 'time_harmonic'
        E1   = sps_getp(Ep, 'E1', 0);
        Edir = sps_unit(sps_getp(Ep, 'Edir', [1;0;0]));
        om   = sps_getp(Ep, 'omega', 0);
        ph0  = sps_getp(Ep, 'phase', 0);
        E    = cfg.E0(:) + Edir*E1*cos(om*t + ph0);

    case 'circular'
        E1  = sps_getp(Ep, 'E1', 0);
        om  = sps_getp(Ep, 'omega', 0);
        ph0 = sps_getp(Ep, 'phase', 0);
        [e1, e2] = sps_perp_basis(b0, b0);
        e1u = sps_getp(Ep, 'e1', []);
        if ~isempty(e1u)
            b0h = sps_unit(b0);
            if norm(b0h) == 0, b0h = [0;0;1]; end
            e1u = e1u(:);
            e1  = sps_unit(e1u - dot(e1u,b0h)*b0h);
            if norm(e1) == 0
                [e1, e2] = sps_perp_basis(b0, b0);
            else
                e2 = cross(b0h, e1);
            end
        end
        E = cfg.E0(:) + E1*(e1*cos(om*t + ph0) + e2*sin(om*t + ph0));

    case 'convective'
        U  = sps_getp(Ep, 'U', [4e5; 0; 0]);  U = U(:);
        Ua = sps_getp(Ep, 'Uamp', 0);          % U 的时变幅度(相对值)
        Uo = sps_getp(Ep, 'Uomega', 0);
        if Ua ~= 0
            U = U * (1 + Ua*sin(Uo*t));
        end
        E = -cross(U, B);

    case 'inductive_wave'
        kv = sps_getp(Bp, 'k', [0;0;0]);  kv = kv(:);
        om = sps_getp(Bp, 'omega', 0);
        kn = norm(kv);
        if kn > 0 && om ~= 0
            E = -(om/kn) * cross(kv/kn, Bwave);
        else
            E = zeros(3,1);
        end

    case 'parallel_wave'
        kv  = sps_getp(Bp, 'k', [0;0;0]);  kv = kv(:);
        om  = sps_getp(Bp, 'omega', 0);
        E1  = sps_getp(Ep, 'E1', 0);
        ph0 = sps_getp(Ep, 'phase', 0);
        ph  = dot(kv, xs) - om*t + ph0;
        E   = E1 * bh * cos(ph);

    case 'turbulence_inductive'
        % 剪切 Alfven 湍流: E = -VA*(b0 x dB),逐模满足 Faraday 定律
        % (omega = VA*k_par),且 E . b0 = 0,不引入平行电场。
        b0h = sps_unit(B0);
        VA  = sps_getp(Bp, 'VA', []);
        if isempty(VA) || norm(b0h) == 0
            warning('sps_fields:noVA', ...
                'turbulence_inductive 需要先经 sps_prepare 生成湍流模与 VA,当前返回 E = 0。');
            E = zeros(3,1);
        else
            E = -VA*cross(b0h, Bwave);
        end

    case 'custom'
        if ~isfield(cfg,'Efun') || isempty(cfg.Efun)
            error('sps_fields:customE', 'Emodel=''custom'' 需要 cfg.Efun = @(t,x) ...');
        end
        E = cfg.Efun(t, x);
        E = E(:);

    otherwise
        error('sps_fields:Emodel', '未知电场模型 "%s"', cfg.Emodel);
end

% ---- E 的整体时间调制 ----
if isfield(cfg,'Emod') && isstruct(cfg.Emod)
    a = sps_getp(cfg.Emod, 'amp', 0);
    if a ~= 0
        E = E * (1 + a*sin(sps_getp(cfg.Emod,'omega',0)*t + sps_getp(cfg.Emod,'phase',0)));
    end
end

aux = struct('Bmag', Bn, 'bhat', bh, 'Bwave', Bwave, 'xs', xs);

end
