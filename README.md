# 空间等离子体环境粒子模拟(Space Plasma Single-Particle Simulator)

一组 MATLAB 函数,用来模拟**空间等离子体环境中的带电粒子运动**。
核心出发点是洛伦兹力:

```
m dv/dt = q (E + v × B)          dx/dt = v
```

- **E、B 都可以随时间变化**(平面波、ULF 扰动、|B| 的整体压缩/呼吸、时谐电场……)
- **E、B 都可以随空间变化**(线性梯度、偶极场、电流片、磁镜、Parker 螺旋、湍流……)
- **场线可以是"直的"也可以是"弯的"**:`uniform` / `linear_gradient` 是直线场线(无曲率漂移),
  `dipole` / `mirror` / `harris` / `turbulent` 是弯曲场线(出现梯度-B 与曲率漂移)。
- 所有初始条件(位置、速率、投掷角、回旋相位、漂移速度、粒子种类、步长……)既可以**从预设/菜单选**,
  也可以**手动输入任意数值**。

---

## 1. 快速开始

把本目录设为 MATLAB 当前目录(或 `addpath` 该目录),然后:

```matlab
space_plasma_sim                       % 【交互菜单】选场景 -> 改参数(可手动输入) -> 运行
space_plasma_sim('list')               % 列出全部预设算例
list = sps_preset_list();              % 取预设清单(编程用)
space_plasma_sim('preset','dipole_belt')   % 直接运行某个预设
space_plasma_sim('demo')               % 多算例演示:预设总览 + 4 组参数扫描
space_plasma_sim('selftest')           % 物理/数值自检(39 项检查,当前全部通过)
```

在脚本里手动构造配置(推荐方式,便于批量研究):

```matlab
cfg = sps_presets('dipole_belt');      % 从预设出发
cfg.init.x0    = [5*6.3712e6; 0; 0];   % 手动改到 L = 5
cfg.init.alpha = 60;                   % 手动改投掷角
cfg.numerics.ngyro = 2000;             % 模拟时长(回旋周期数)
out = sps_run(cfg);                    % 运行,自动出诊断报告 + 图
```

或者从零构造(所有字段及其含义见 `sps_default_config.m`):

```matlab
cfg = sps_default_config();
cfg.Bmodel = 'harris';
cfg.Bparams = struct('Bx',20e-9,'Bn',3e-9,'L',1e6);
cfg.Emodel = 'uniform'; cfg.E0 = [1e-3;0;0];
cfg.init.mode = 'gc';                   % 用 速率/投掷角/回旋相位
cfg.init.x0 = [0;0;0]; cfg.init.speed = 5e5; cfg.init.alpha = 60;
out = sps_run(cfg);
```

---

## 2. 文件说明

| 文件 | 作用 |
|---|---|
| `space_plasma_sim.m` | 主入口(菜单 / 预设 / 演示 / 自检) |
| `sps_default_config.m` | **全部可选项与默认值**都在这里,先读它 |
| `sps_presets.m` | 10 个真实空间环境的预设算例 |
| `sps_preset_list.m` | 预设清单(名称/场景/物理/场结构) |
| `sps_config_menu.m` | 文本交互菜单,任何数值都能手动输入 |
| `sps_fields.m` | **场模型库**:E(t,x) 与 B(t,x) 的所有时空变化形式 |
| `sps_boris.m` | Boris 推进器(可选相对论版本) |
| `sps_run.m` | 主求解器:时间积分 + 存储 |
| `sps_analysis.m` | 诊断:守恒量、导心漂移、回旋/反弹周期 |
| `sps_drifts.m` | 导心漂移各解析分量(E×B、∇B、曲率、极化、平行) |
| `sps_plots.m` | 轨道/速度/投掷角/守恒量/漂移对比图 |
| `sps_prepare.m` | 校验配置、自动定时间步长、自适应子步设置 |
| `sps_initial_state.m` | 把 速率/投掷角/回旋相位 化成三维初速度 |
| `sps_make_turbulence.m` | 生成多模 Alfvén 湍流场(不可压、`ω=V_A k_∥`、可复现) |
| `sps_demo.m` | 预设总览 + 参数扫描演示 |
| `sps_selftest.m` | 自检(解析对比、div B、Faraday 一致性、共振散射) |
| `sps_const.m` / `sps_species.m` / `sps_unit.m` / `sps_getp.m` / `sps_v_from_E.m` / `sps_perp_basis.m` / `sps_print_config.m` | 常数、粒子、向量工具 |

单位:全部 SI(位置 m、速度 m/s、时间 s、磁场 T、电场 V/m);投掷角/相位用**度**。

---

## 3. 磁场模型(`cfg.Bmodel` 与 `cfg.Bparams`)

| `Bmodel` | 场位形 | 场线弯曲? | 主要参数 |
|---|---|---|---|
| `uniform` | 均匀恒定 `B = B0` | 直 | `cfg.B0` |
| `linear_gradient` | `|B| = Bmag(1+amp·(gdir·x)/L)`,方向 `Bdir` | 直(平行场线) | `Bdir,Bmag,grad_dir,grad_amp,L,origin` |
| `dipole` | 偶极场 `B = μ0/4π·(3(m·r̂)r̂−m)/r³` | **弯** | `moment`(或 `eq_field`+`R_planet`)、`origin` |
| `harris` | 电流片 `B = Bx·tanh(s)·t̂ + Bg·ĝ + Bn·n̂` | **强弯(可反平行)** | `Bx,Bg,Bn,L,nhat,tangent` |
| `mirror` | 有限长磁镜 `B_z = Bmag[1+(R−1)sin²(πz/2L)]`,`B_ρ=−(ρ/2)∂B_z/∂z` | **会聚** | `Bmag,L,ratio,axis` |
| `wave` | 背景场 + 平面波(圆/线偏振) | 直 | `k,omega,B1,polarization('RH'/'LH'/'linear'),phase0,e1` |
| `parker` | Parker 螺旋太阳风场 | 弯(螺旋) | `r0,B_r0,Vsw,OmegaSun,axis` |
| `turbulent` | 背景场 + N 个 Alfvén 模(不可压,随机相位,可复现),`ω = V_A·k_∥` | 随机弯 | `Nmodes,k_min,k_max,slope,dB_rms,n0,seed` |
| `custom` | 用户函数 `cfg.Bfun = @(t,x) [Bx;By;Bz]` | 任意 | — |

**时间变化**另有两种方式,可与任意空间模型叠加:

1. 波/湍流自带 `omega`(相位随时间演化);
2. 全局幅值调制 `cfg.Bmod = struct('amp',a,'omega',w,'phase',p)`:
   `B(t) = B_model(t,x) · (1 + a·sin(w t + p))`(模拟磁压缩/呼吸模)。

> ⚠ `cfg.Bmod` 是**运动学**调制:它只改变 |B|,不附带 Faraday 定律要求的感应电场,
> 因此单独使用时不是自洽的电磁场(此时磁矩会按 `μ ∝ 1/|B|` 变化)。
> 需要自洽的时变场时,用 `Bmodel='wave'` + `Emodel='inductive_wave'`
> (演示 `sps_demo` 的第 4 组扫描就是这个用法:ω/Ω_c 从 0.05 到 3,
> 可以看到磁矩从守恒(绝热)过渡到不守恒(非绝热))。

> 关于 `div B = 0`:除 `custom` 外所有模型都保证无散度;
> `linear_gradient` 会自动把 `grad_dir` 投影到垂直 `Bdir` 的方向(平行梯度无法自洽)。
> 自检脚本对每个模型都用中心差分验证了这一点。

---

## 4. 电场模型(`cfg.Emodel` 与 `cfg.Eparams`)

| `Emodel` | 表达式 | 参数 |
|---|---|---|
| `zero` | `E = 0` | — |
| `uniform` | `E = E0` | `cfg.E0` |
| `time_harmonic` | `E = E0 + E1·cos(ωt+φ)·Edir` | `E1,Edir,omega,phase` |
| `circular` | 垂直 B 的旋转电场 | `E1,omega,phase,e1` |
| `convective` | `E = −U × B`(太阳风/磁层对流) | `U`,`Uamp`,`Uomega` |
| `inductive_wave` | Faraday 定律 `E = −(ω/|k|) k̂ × δB` | 配合 `Bmodel='wave'` |
| `parallel_wave` | `E = E1·b̂·cos(k·x−ωt)`(平行电场) | `E1`,`phase` |
| `turbulence_inductive` | 各湍流模的感应电场 `E = −V_A × δB`(理想 Alfvén 波,**无平行电场**) | 配合 `Bmodel='turbulent'` |
| `custom` | `cfg.Efun = @(t,x) [Ex;Ey;Ez]` | — |

同样可叠加全局幅值调制 `cfg.Emod`。

---

## 5. 初始条件(`cfg.init`)

| 字段 | 含义 |
|---|---|
| `mode = 'gc'` | 用 速率 + 投掷角 + 回旋相位 给出速度(推荐,**相对局地磁场**定义) |
| `mode = 'cart'` | 直接给三维速度矢量 |
| `x0` | 初始位置 [m] |
| `speed` | 速率 [m/s] |
| `alpha` | 投掷角 [deg];`v_∥ = v·cosα`,`v_⊥ = v·sinα` |
| `gyrophase` | 回旋相位 [deg] |
| `vpar`,`vperp` | 若给有限值则直接使用(忽略 `speed/alpha`) |
| `vdrift` | 附加的整体漂移速度 [m/s] |
| `v0` | `mode='cart'` 时的初速度 |

粒子:

| `cfg.species` | 说明 |
|---|---|
| `proton` / `electron` / `alpha` / `O+` | 内置质量电荷 |
| `custom` | 需自己给 `cfg.q`(带符号,C)与 `cfg.m`(kg) |

`cfg.relativistic = true` 时使用相对论 Boris 推进器(动量 `u = γv`),
适合辐射带高能电子等(预设 `whistler_electron` 就是 53 keV、0.42c 的电子)。

---

## 6. 十种情况的初始条件(预设算例)

> 下表中的 `Ω_c`、`T_c`、`r_g` 是按初始点局地场算的近似值,
> 程序运行时会打印**精确**的配置摘要与诊断报告。

| 预设 | 磁场 | 电场 | 初始位置 | 初速度 | 回旋尺度 | 观察什么 |
|---|---|---|---|---|---|---|
| `gyro_basic` | 均匀 `B=50 nT ẑ` | `E=1 mV/m x̂` | 原点 | 质子 100 km/s, α=90°, φ=0 | `T_c=1.31 s`, `r_g=20.9 km` | 回旋 + `E×B=20 km/s(−ŷ)` |
| `gradient_drift` | `B_z=50 nT(1+0.5x/L)`, `L=1000 km` | 0 | 原点 | 质子 100 km/s, α=90° | `T_c=1.31 s`, `r_g=20.9 km` | ∇B 漂移 ≈0.52 km/s(+ŷ) |
| `dipole_belt` | 地球偶极 `m=8.07e22 A·m²` | 0 | `x=4 R_E` 赤道 | 1 MeV 质子 1.384e7 m/s, α=90° | `B=488 nT`, `T_c=0.135 s`, `r_g=296 km` | **弯曲场线**引起的梯度+曲率漂移(西向) |
| `convection` | 地球偶极 | 晨昏 0.5 mV/m + 100 s ULF 0.2 mV/m(**时变 E**) | `x=6 R_E` 赤道 | 1 keV 质子 4.38e5 m/s, α=90° | `B=144 nT`, `T_c=0.45 s`, `r_g=32 km` | `E×B≈3.5 km/s` 向日内输运 + 极化漂移 |
| `solar_wind` | Parker 螺旋:`r0=10 R_sun`, `B_r0=1.386 µT` → 1 AU 处 `B_r=3 nT, B_φ=−3.06 nT` | 对流 `E=−U×B`,`U=400 km/s` | `x=1 AU` | 质子热速度 30 km/s(α=45°)+ 随太阳风整体流动 400 km/s | `B=4.29 nT`, `T_c=15.3 s`, `r_g=1048 km` | 对流电场使导心随太阳风运动,`E×B` 漂移自动给出流速 |
| `tail_current_sheet` | Harris:`Bx=20 nT`, `L=1000 km`, `Bn=3 nT` | 0 | 电流片中心 | 质子 500 km/s, α=60° | 中心 `B=3 nT`, `T_c=22 s` | Speiser 型混沌轨道、反复穿越中性片 |
| `mirror_trap` | 磁镜:`B0=100 nT`, `R=4`, `L=1000 km` | 0 | 原点(镜心) | 质子 200 km/s, α=45°, φ=45° | `T_c=0.656 s` | 损失锥角 30°:α=45° 被捕获来回反弹 |
| `whistler_electron` | `B0=200 nT` + 哨声波 `f=2 kHz`, `k=−2.51e-4 ẑ`, `δB=10 nT` | 波感应电场(`|E|≈0.5 V/m`) | 原点 | **相对论电子** `v_∥=v_⊥=9e7 m/s`(53 keV) | `T_c=178 µs` | 回旋共振 `ω−k v_∥=Ω_ce` → 投掷角散射 |
| `emic_proton` | `B0=100 nT` + EMIC 波 `f=0.5 Hz`, `k=4.55e-6 ẑ`, `δB=5 nT` | 波感应电场(`|E|≈3.5 mV/m`) | 原点 | 质子 3.23e6 m/s, α=30°(54 keV) | `T_c=0.656 s` | 质子回旋共振 `ω−k v_∥=−Ω_cp` |
| `turbulence` | `B0=5 nT` + 48 个 Alfvén 模,`δB_rms=1.5 nT`,波长 2e6~1e8 m | `E=−V_A×δB`(无平行电场) | 原点 | 质子 100 km/s, α=60° | `T_c=13.1 s`, `r_g=209 km` | 随机场中的扩散、投掷角随机化 |

想快速看"不同情况"的区别:

```matlab
space_plasma_sim('demo')
```

会跑:

1. 全部预设的对比表(`|B|`、`T_c`、`r_g`、回旋平均后的磁矩守恒性、实测 vs 解析漂移);
2. **投掷角扫描**(磁镜:10/25/45/70 度 → 逃逸 vs 捕获、反弹周期);
3. **电场扫描**(0/0.5/1/2 mV/m → `E×B` 漂移速度线性变化);
4. **弯曲 vs 直场**(偶极场 vs 同强度均匀场 → 有无曲率/梯度漂移);
5. **场时变快慢扫描**(ω/Ω_c = 0.05 / 0.5 / 3 → 绝热时磁矩守恒,非绝热时投掷角被散射)。

---

## 7. 数值方法

* **Boris 推进器**:电场半步加速 + 磁场纯旋转,磁场步严格保速度模长,
  长时间积分能量/磁矩守恒性好(等离子体 PIC 的标准做法);可选相对论版本。
* **时间步长自动选取**:`dt = min(T_c, T_wave)/nper`,保证每个回旋周期至少 `nper` 步
  (默认 120),同时分辨最快的场时标(波、调制、湍流模)。
  也可用 `cfg.numerics.dt` / `cfg.numerics.tmax` 直接指定。
* **自适应子步**:粒子走进比参考场强得多的区域时(局地场超过参考场
  `cfg.numerics.substeps_trigger` 倍,默认 2 倍),自动把一步细分成多个子步作为安全网。
  阈值故意设得不敏感:若一有风吹草动就改变子步数,积分格式会随粒子位置变得不对称,
  从而引入虚假的磁矩漂移(数值实验:阈值取 1 时偶极场 `|Δμ|/μ` 会从 2e-3 恶化到 7e-2)。
* **逃逸/失稳保护**:`|B|` 超阈值或 `|x|` 超界时提前终止并给出原因,
  出现 NaN 也会被捕获(预设里都能看到 `out.escaped`)。
* **诊断**:磁矩 `μ = p_⊥²/(2mB)`、动能、投掷角、瞬时导心位置
  `R_gc = x + m(v×B)/(qB²)`,并给出六个漂移分量的解析值:

```
v_par·b̂                          平行运动
E×B/B²                           电场漂移
μ(B×∇B)/(qB³)                    梯度-B 漂移
m v_∥²(κ×B)/(qB²),  κ=(b̂·∇)b̂    曲率漂移(空间差分求 κ)
m/(qB²)·dE_⊥/dt                  极化漂移
```

  注意:粒子位置上直接算出的 `μ` 会有一个 `O(r_g|∇B|/B)` 量级的回旋摆动
  (偶极场 L=4 的 1 MeV 质子约 ±6.6%,理论上 `2ε`,ε = r_g|∇B|/B = 3.5%),
  这是真实物理而非数值误差。程序同时给出**回旋平均后**的 `μ`,它才是真正的
  第一绝热不变量。经数值实验检验(偶极场 L=4、600 个回旋周期):
  回旋平均后的 `|Δμ|/μ ≈ 2×10⁻³`,并且随步长减小不再变化(已收敛),
  残余量按 `ε²` 标度,属于有限拉莫半径修正。

---

## 8. 自检

```matlab
space_plasma_sim('selftest')
```

包含(全部与解析结果对比,逐项 PASS/FAIL,共 39 项,当前版本全部通过):

1. 均匀场纯回旋:动能/磁矩守恒、回旋周期、回旋半径;
2. `E×B` 漂移的大小与方向;
3. 偶极场:回旋平均磁矩的绝热守恒、质子西向漂移、漂移速度量级;
4. 磁镜:损失锥内逃逸 / 损失锥外捕获、反弹周期;
5. 8 个场模型的 `div B ≈ 0`(中心差分数值检验);
6. 感应电场量级 `|E| = (ω/k)·δB`;
7. **Faraday 定律一致性** `∇×E + ∂B/∂t ≈ 0`(波场与湍流场),以及理想 Alfvén 波
   的**无平行电场** `E·b̂₀ = 0`;
8. 相对论推进:纯磁场下 `γ` 守恒;
9. 10 个预设都能正常运行;
10. 波-粒回旋共振使投掷角明显变化(哨声波-电子、EMIC-质子)。

---

## 9. 适用范围与限制

* **试验粒子(test-particle)模型**:只跟踪单个粒子,不求解等离子体自洽响应
  (没有 PIC 的场-粒子反馈、没有集体效应)。
* **无碰撞**:未包含库仑碰撞/波-粒扩散的统计处理;波-粒相互作用是"单粒子在给定波场中运动"。
* 场模型是理想化解析位形(例如偶极场是纯偶极,没有地磁尾拉伸、没有屏蔽效应)。
* 强场区、零场点(`harris` 中心且 `Bn=0`、偶极奇点)附近导心近似失效,
  程序会给出提示或提前终止。
* 需要 MATLAB R2016b 及以上(使用了函数文件内的局部函数);推荐 R2020a+ 以正确显示中文注释。

## 10. 扩展方式

* **加新场模型**:在 `sps_fields.m` 的 `switch` 里加一个 `case`,并根据需要补充 `cfg.Bparams` 字段;
  菜单项在 `sps_config_menu.m` 的 `edit_B` 中加一行即可。
* **完全自定义**:`cfg.Bmodel='custom'; cfg.Bfun=@(t,x)[...];`
  `cfg.Emodel='custom'; cfg.Efun=@(t,x)[...];`
* **批量扫描**:写个 `for` 循环改 `cfg.init.*` / `cfg.*0` / `cfg.Bparams.*`,
  用 `cfg.verbose=false; cfg.plot=false;` 静默运行,再从 `out.summary` 里取量做统计。

---

## 11. 许可证

MIT License,见 [LICENSE](LICENSE)。可以自由使用、修改、再分发,保留版权声明即可。

## 12. 引用 / 致谢

如果这套代码对你有帮助,欢迎在论文或报告中引用本仓库。
物理模型与数值方法参考:

* Boris, J. P. (1970). *Relativistic plasma simulation—optimization of a hybrid code.*
  Proc. 4th Conf. Numerical Simulation of Plasmas.
* Northrop, T. G. (1963). *The Adiabatic Motion of Charged Particles.* Interscience.
* Parker, E. N. (1958). *Dynamics of the Interplanetary Gas and Magnetic Fields.*
  ApJ, 128, 664.
* Baumjohann, W. & Treumann, R. A. (1996). *Basic Space Plasma Physics.* Imperial College Press.
