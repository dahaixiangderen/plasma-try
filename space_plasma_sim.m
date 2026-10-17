function out = space_plasma_sim(action, arg)
%SPACE_PLASMA_SIM  空间等离子体环境粒子模拟 —— 主入口
%
%   在 MATLAB 中把本文件所在目录设为当前目录(或将目录加入路径),然后:
%
%     space_plasma_sim                     % 【交互菜单】选场景/改参数/手动输入数值/运行
%     space_plasma_sim('list')             % 列出全部预设算例
%     space_plasma_sim('preset','dipole_belt')   % 直接运行某个预设
%     space_plasma_sim('demo')             % 多算例演示:预设对比 + 参数扫描
%     space_plasma_sim('selftest')         % 运行物理/数值自检
%     cfg = space_plasma_sim('config')     % 只打开菜单拿配置,不运行
%     out = space_plasma_sim(cfg)          % 直接运行自己构造的 cfg
%
%   输出 out 见 sps_run 的说明(轨道、场、守恒量、诊断汇总)。
%
%   物理模型与全部选项见 README.md;场模型实现见 sps_fields.m。

here = fileparts(mfilename('fullpath'));
if ~isempty(here)
    addpath(here);
end

if nargin < 1 || isempty(action)
    action = 'menu';
end

% 直接给配置结构体
if isstruct(action)
    out = sps_run(action);
    return;
end

switch lower(strtrim(action))

    case 'menu'
        out = [];
        while true
            cfg = sps_config_menu();
            if isempty(cfg)
                fprintf('已取消。\n');
                return;
            end
            out = sps_run(cfg);
            fprintf('\n本次运行完成(算例: %s)。\n', cfg.tag);
            s = input('是否继续调整并再算一次? [y/N]: ', 's');
            if isempty(strtrim(s)) || ~any(strcmpi(strtrim(s), {'y','yes','1','是'}))
                return;
            end
        end

    case 'preset'
        if nargin < 2 || isempty(arg)
            error('space_plasma_sim:needname', '用法: space_plasma_sim(''preset'',''dipole_belt'')');
        end
        cfg = sps_presets(arg);
        out = sps_run(cfg);

    case 'list'
        list = sps_preset_list();
        fprintf('\n================ 预设算例 ================\n');
        fprintf('%-20s %-32s %-24s %s\n', '名称', '场景', '主要物理', '场结构');
        for i = 1:numel(list)
            fprintf('%-20s %-32s %-24s %s\n', ...
                    list(i).name, list(i).title, list(i).physics, list(i).field);
        end
        fprintf('==========================================\n');
        fprintf('运行方式: space_plasma_sim(''preset'',''<名称>'')\n');
        out = list;

    case 'demo'
        sps_demo();
        out = [];

    case 'selftest'
        ok = sps_selftest();
        out = ok;

    case 'config'
        out = sps_config_menu();

    otherwise
        % 也允许直接传预设名: space_plasma_sim('dipole_belt')
        try
            cfg = sps_presets(action);
        catch
            error('space_plasma_sim:usage', ...
                  ['未知命令 "%s"。可用: menu / list / preset / demo / selftest / config,', ...
                   '或直接给预设名。'], action);
        end
        out = sps_run(cfg);
end

end
