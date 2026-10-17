function c = sps_const()
%SPS_CONST  物理常数与常用空间尺度(全部为国际单位制 SI)
%
%   c = SPS_CONST() 返回结构体 c,字段含义见下:
%     e,mp,me,malpha,mO : 电荷与各粒子质量
%     mu0,eps0,c_light  : 真空磁导率/介电常数/光速
%     kB                : 玻尔兹曼常数
%     Re, AU            : 地球半径, 天文单位
%     OmegaSun          : 太阳自转角速度
%     B_eq_earth        : 地球赤道表面磁场强度

c.e          = 1.602176634e-19;   % 基本电荷 [C]
c.mp         = 1.67262192369e-27; % 质子质量 [kg]
c.me         = 9.1093837015e-31;  % 电子质量 [kg]
c.malpha     = 6.6446573357e-27;  % alpha 粒子质量 [kg]
c.mO         = 2.6561e-26;        % 氧离子 O+ 质量 [kg]

c.mu0        = 4*pi*1e-7;         % 真空磁导率 [H/m]
c.eps0       = 8.8541878128e-12;  % 真空介电常数 [F/m]
c.c_light    = 299792458;         % 光速 [m/s]
c.kB         = 1.380649e-23;      % 玻尔兹曼常数 [J/K]

c.Re         = 6.3712e6;          % 地球平均半径 [m]
c.AU         = 1.495978707e11;    % 天文单位 [m]
c.Rsun       = 6.957e8;           % 太阳半径 [m]
c.OmegaSun   = 2.865e-6;          % 太阳自转角速度 [rad/s]

c.B_eq_earth = 3.12e-5;           % 地球赤道表面磁场 [T] (~0.312 G)
c.m_earth    = 8.07e22;           % 地球偶极矩大小 [A*m^2]

end
