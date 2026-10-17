function [q, m, label] = sps_species(cfg)
%SPS_SPECIES  根据 cfg.species 返回电荷 q [C]、质量 m [kg] 和中文名称
%
%   支持: 'proton' (H+), 'electron' (e-), 'alpha' (He2+), 'O+' (氧离子),
%         'custom' (此时必须给 cfg.q 与 cfg.m)
%
%   注意: q 带符号(电子为负),回旋频率 Omega = qB/m 也随之带符号。

c = sps_const();

switch lower(strtrim(cfg.species))
    case {'proton','p','h+','ion','hydrogen'}
        q = c.e;      m = c.mp;     label = '质子 H+';
    case {'electron','e','e-','electron-'}
        q = -c.e;     m = c.me;     label = '电子 e-';
    case {'alpha','he2+','he++','helium'}
        q = 2*c.e;    m = c.malpha; label = 'alpha 粒子 He2+';
    case {'o+','oxygen','o'}
        q = c.e;      m = c.mO;     label = '氧离子 O+';
    case {'custom','user'}
        q = cfg.q;    m = cfg.m;    label = '自定义粒子';
    otherwise
        error('sps_species:unknown', ...
              '未知粒子种类 "%s" (可选: proton/electron/alpha/O+/custom)', cfg.species);
end

if isempty(q) || isempty(m) || ~isfinite(q) || ~isfinite(m) || m <= 0
    error('sps_species:badcustom', ...
          'species=''custom'' 时必须提供有限的 cfg.q(电荷) 与正的 cfg.m(质量)。');
end

end
