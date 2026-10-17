function [x0, v0] = sps_initial_state(cfg)
%SPS_INITIAL_STATE  由配置解析出初始位置 x0 [m] 与速度 v0 [m/s]
%
%   init.mode = 'cart' : 直接使用 init.x0 与 init.v0
%   init.mode = 'gc'   : 由 速率/投掷角/回旋相位 构造速度
%        v_par  = v cos(alpha)   (alpha 以局部磁场方向为参考)
%        v_perp = v sin(alpha)
%        v0     = v_par*bhat + v_perp*(cos(phi)*e1 + sin(phi)*e2) + vdrift
%      其中 (e1,e2,bhat) 为右手垂直基,phi 为回旋相位。
%      若 init.vpar / init.vperp 为有限值,则直接采用它们。

init = cfg.init;
x0   = init.x0(:);
init.vdrift = init.vdrift(:);

switch lower(init.mode)

    case 'cart'
        v0 = init.v0(:);

    case 'gc'
        [~, B] = sps_fields(0, x0, cfg);
        bh = sps_unit(B);
        if norm(bh) == 0
            bh = [0; 0; 1];                % 零场退化为以 z 为参考
        end
        [e1, e2] = sps_perp_basis(bh, bh);

        vpar  = init.vpar;
        vperp = init.vperp;
        sp    = sps_getp(init, 'speed', 1e5);
        al    = sps_getp(init, 'alpha', 90)*pi/180;

        if ~isfinite(vpar),  vpar  = sp*cos(al); end
        if ~isfinite(vperp), vperp = sp*sin(al); end

        ph = sps_getp(init, 'gyrophase', 0)*pi/180;
        v0 = vpar*bh + vperp*(cos(ph)*e1 + sin(ph)*e2) + init.vdrift;

    otherwise
        error('sps_initial_state:mode', ...
              '未知 init.mode "%s" (可选 ''gc'' 或 ''cart'')', init.mode);
end

end
