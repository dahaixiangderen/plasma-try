function v = sps_v_from_E(E_eV, m)
%SPS_V_FROM_E  由动能(非相对论)换算速率: v = sqrt(2E/m)
%
%   E_eV : 动能 [eV] (标量或数组)
%   m    : 粒子质量 [kg]
%   v    : 速率 [m/s]

c = sps_const();
v = sqrt(2 * E_eV * c.e / m);

end
