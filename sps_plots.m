function sps_plots(out)
%SPS_PLOTS  绘制模拟结果:轨道、速度与投掷角、守恒量与漂移对比
%
%   SPS_PLOTS(out)  其中 out 为 sps_run 的返回值

cfg = out.cfg;
nt  = numel(out.t);
if nt < 2
    warning('sps_plots:tooshort', '数据点太少,跳过绘图。');
    return;
end

% 抽样,避免绘图点过多
step = max(1, floor(nt/4000));
idx  = 1:step:nt;
if idx(end) ~= nt, idx(end+1) = nt; end

t  = out.t(idx);
x  = out.x(:, idx);
v  = out.v(:, idx);
gc = out.gc(:, idx);
Bm = out.Bmag(idx);
Em = out.Emag(idx);
pa = out.pitch(idx);

% 轨道上的平行/垂直速度(直接用已存场,避免重复求场)
Bd  = out.B(:, idx);
bn  = sqrt(sum(Bd.^2, 1));
bh  = Bd ./ max(bn, realmin);
vparv  = sum(v.*bh, 1);
vperp2 = max(sum(v.^2, 1) - vparv.^2, 0);
vperpv = sqrt(vperp2);

% 长度单位
span = max(out.x, [], 2) - min(out.x, [], 2);
[Lu, LuName] = length_unit(max(span));

% 初值处的垂直基(用于回旋面投影)
bh1 = sps_unit(out.B(:,1));
[e1, e2] = sps_perp_basis(bh1, bh1);

% ======================================================================
% 图 1: 轨道
% ======================================================================
figure('Name', sprintf('轨道 - %s', cfg.tag), 'Color', 'w');
set(gcf, 'Units','normalized', 'OuterPosition', [0.02 0.06 0.54 0.86]);

subplot(2,2,[1 3]);
hold on; grid on; box on;
hFL = [];
if ~strcmpi(cfg.Bmodel,'custom')
    hFL = draw_fieldlines(out, cfg, Lu);
end
hSc = scatter3(x(1,:)/Lu, x(2,:)/Lu, x(3,:)/Lu, 6, Bm*1e9, 'filled');
hGC = plot3(gc(1,:)/Lu, gc(2,:)/Lu, gc(3,:)/Lu, 'k-', 'LineWidth', 0.8);
hS  = plot3(x(1,1)/Lu, x(2,1)/Lu, x(3,1)/Lu, 'go', 'MarkerFaceColor','g', 'MarkerSize',7);
hE  = plot3(x(1,end)/Lu, x(2,end)/Lu, x(3,end)/Lu, 'rs', 'MarkerFaceColor','r', 'MarkerSize',7);
cb = colorbar; cb.Label.String = '|B| (nT)';
xlabel(sprintf('x (%s)', LuName)); ylabel(sprintf('y (%s)', LuName));
zlabel(sprintf('z (%s)', LuName));
title(sprintf('%s  轨道(颜色=|B|)', cfg.tag));
if isempty(hFL)
    legend([hSc hGC hS hE], {'粒子(真实)','导心','起点','终点'}, 'Location','best');
else
    legend([hFL hSc hGC hS hE], {'磁力线','粒子(真实)','导心','起点','终点'}, 'Location','best');
end
view(35, 20); axis equal;

subplot(2,2,2);
hold on; grid on; box on;
rm = mean(x, 2);
plot((e1.'*(x - rm))/Lu, (e2.'*(x - rm))/Lu, 'b-');
plot((e1.'*(x(:,1)-rm))/Lu, (e2.'*(x(:,1)-rm))/Lu, 'go', 'MarkerFaceColor','g');
if Lu >= 1e3
    title(sprintf('垂直磁场平面上的投影(回旋轨道, %s)', LuName));
else
    title('垂直磁场平面上的投影(回旋轨道)');
end
xlabel(sprintf('e_1 (%s)', LuName)); ylabel(sprintf('e_2 (%s)', LuName));
axis equal;

subplot(2,2,4);
hold on; grid on; box on;
plot(x(1,:)/Lu, x(2,:)/Lu, 'b-');
plot(gc(1,:)/Lu, gc(2,:)/Lu, 'r-', 'LineWidth', 1.2);
xlabel(sprintf('x (%s)', LuName)); ylabel(sprintf('y (%s)', LuName));
title('x-y 平面(可见漂移)'); legend({'粒子','导心'},'Location','best');
axis equal;

% ======================================================================
% 图 2: 速度、投掷角、场强
% ======================================================================
figure('Name', sprintf('速度与投掷角 - %s', cfg.tag), 'Color','w');
set(gcf, 'Units','normalized', 'OuterPosition', [0.58 0.52 0.40 0.42]);

subplot(3,1,1);
plot(t, vparv/1e3, 'r-', t, vperpv/1e3, 'b-', t, sqrt(sum(v.^2,1))/1e3, 'k--');
grid on; legend({'v_\parallel','v_\perp','|v|'},'Location','best');
ylabel('速度 (km/s)'); title('速度分量');

subplot(3,1,2);
plot(t, pa, 'm-'); grid on;
ylabel('投掷角 \alpha (deg)'); ylim([0 180]);
title(sprintf('投掷角随时间变化(范围 %.2f^o ~ %.2f^o)', min(pa), max(pa)));

subplot(3,1,3);
semilogy(t, max(Bm*1e9, realmin), 'b-', t, max(Em*1e3, realmin), 'r-'); grid on;
legend({'|B| (nT)','|E| (mV/m)'},'Location','best');
xlabel('时间 (s)'); title('轨道上的场强');

% ======================================================================
% 图 3: 守恒量与漂移
% ======================================================================
figure('Name', sprintf('守恒量与漂移 - %s', cfg.tag), 'Color','w');
set(gcf, 'Units','normalized', 'OuterPosition', [0.58 0.05 0.40 0.42]);

subplot(2,2,1);
plot(out.t(idx), (out.KE(idx) - out.KE(1))/max(abs(out.KE(1)),realmin), 'b-');
grid on; xlabel('t (s)'); ylabel('\Delta KE / KE_0');
title('动能(相对)变化');

subplot(2,2,2);
plot(out.t(idx), (out.mu(idx) - out.mu(1))/max(abs(out.mu(1)),realmin), 'r-');
grid on; xlabel('t (s)'); ylabel('\Delta \mu / \mu_0');
title('第一绝热不变量 \mu 的漂移');

subplot(2,2,[3 4]);
hold on; grid on; box on;
s = out.summary;
if isfield(s.drift_initial,'vExB')
    d = s.drift_initial;
    vals = [norm(d.vExB), norm(d.vgradB), norm(d.vcurv), norm(d.vparB), norm(d.vpol), norm(d.vtotal)];
    hb = bar(max(vals, 1e-12));
    set(gca, 'YScale','log');
    set(gca, 'XTickLabel', {'E\timesB','\nablaB','曲率','平行','极化','合成'});
    if isfinite(s.drift_speed_measured) && s.drift_speed_measured > 0
        hy = yline(s.drift_speed_measured, 'k--', 'LineWidth', 1.2);
        legend([hb hy], {'解析分量','实测 |v_{gc}|'}, 'Location','best');
        title(sprintf('导心漂移分量(合成 %.3g m/s,实测 %.3g m/s)', ...
              norm(d.vtotal), s.drift_speed_measured));
    else
        legend(hb, {'解析分量'}, 'Location','best');
        title('导心漂移分量(初始点解析值)');
    end
end
ylabel('|v| (m/s)');

end

% ======================================================================
function [Lu, name] = length_unit(span)
%LENGTH_UNIT  根据轨迹尺度选择合适的长度单位

c = sps_const();
if span >= 0.05*c.AU
    Lu = c.AU;   name = 'AU';
elseif span >= 2*c.Re
    Lu = c.Re;   name = 'R_E';
elseif span >= 1e3
    Lu = 1e3;    name = 'km';
else
    Lu = 1;      name = 'm';
end

end

% ======================================================================
function hFirst = draw_fieldlines(out, cfg, Lu)
%DRAW_FIELDLINES  在轨道附近画几条磁力线(直观展示场的弯曲程度)

x0   = out.x(:,1);
span = max(max(out.x,[],2) - min(out.x,[],2));
rg   = out.cfg.scales.rg0;
d    = max([span/6, 4*rg, 1e3]);

[~, B0] = sps_fields(out.t(1), x0, cfg);
b0 = sps_unit(B0);
if norm(b0) == 0
    hFirst = [];
    return;
end
[e1, e2] = sps_perp_basis(b0, b0);

starts = {x0, x0 + 2*d*e1, x0 - 2*d*e1, x0 + 2*d*e2, x0 - 2*d*e2};
Ltot = max(6*d, 8*rg);
np   = 150;

hFirst = [];
for i = 1:numel(starts)
    for sgn = [1 -1]
        P = trace_fline(starts{i}, cfg, sgn*Ltot, np);
        if isempty(P)
            continue;
        end
        h = plot3(P(1,:)/Lu, P(2,:)/Lu, P(3,:)/Lu, '-', ...
                  'Color', [0.55 0.55 0.55], 'LineWidth', 0.5);
        if isempty(hFirst)
            hFirst = h;
        end
    end
end

end

% ======================================================================
function P = trace_fline(xs, cfg, Ltot, npts)
%TRACE_FLINE  沿磁场方向按弧长积分出一条磁力线(RK2)

x  = xs(:);
P  = zeros(3, npts);
ds = Ltot/npts;
for i = 1:npts
    [~, B1] = sps_fields(0, x, cfg);
    b1 = sps_unit(B1);
    if norm(b1) == 0
        P = P(:, 1:max(i-1,1));
        return;
    end
    [~, B2] = sps_fields(0, x + 0.5*ds*b1, cfg);
    b2 = sps_unit(B2);
    if norm(b2) == 0, b2 = b1; end
    x = x + ds*b2;
    if ~all(isfinite(x))
        P = P(:, 1:max(i-1,1));
        return;
    end
    P(:,i) = x;
end

end
