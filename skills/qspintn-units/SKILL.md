---
name: qspintn-units
description: Convert QSpinTN Cm, susceptibility and M-H data between natural, SI and CGS units; identify per-spin, molar and mass normalization from solver outputs and experimental metadata.
---

# Cm、磁化率与 M–H 的单位换算

用于数据初始化、绘图、模型输出对比和结果报告。以下源码路径均相对 **QSpinTN.ml 根目录**。先确定物理量、输出层级、场方向及归一化，再换算数值；不能仅凭变量名 `M` 或标签 `cm³/mol` 判断单位。

## 输出层级与适用范围

机器可读契约仍为 [thermal_units.json](../contracts/thermal_units.json)，标识 `qspintn_chi_si_cm3_v1`。它描述 `ManyBodySolver/GetResult.m` 的物理热力学输出；本 skill 对 M–H 的说明不扩展该 JSON 的结构，也不表示高层接口已支持 M–H 拟合。

| 来源 | 数值含义与需要的处理 |
| --- | --- |
| `UtilityFunc/GetUnitCon.m` | `Model.ES > 0`，单位 K；`T = ES/beta`；`Field.B` 以 T 输入 |
| `ManyBodySolver/ED/ED_func/ED_Cm.m` | 返回全系统自然单位热容 `C`；转每摩尔磁性自旋需 `C/L * R` |
| `ManyBodySolver/ED/ED_func/ED_chi.m` | `M` 是所传算符的热平均；总磁矩算符所得结果需除以 `L` 才为每自旋值 |
| `ManyBodySolver/ED/GetEDRslt.m` | 已将 `Cm` 和 `M` 除以 `Para.L`；不要再次除以 `L` |
| `ManyBodySolver/GetResult.m` 的 `Cm_l` | `real(Rslt.Cm) * R`，J/(mol K)，每摩尔磁性自旋 |
| 同文件的 `Chi_l` | **SI 约定的 cm³/mol**，每摩尔磁性自旋；不是通常实验 CGS cm³/mol |
| 同文件 Chi 分支的 `Rslt_exp_unit.M` | 复制已在 log(beta) 上经 spline 重采样的 `real(Rslt.M)`；在含 g 的每自旋磁矩约定下是 μB/spin 的数值，**未乘 SI/CGS 磁矩换算因子** |

直接调用其他 solver 时核对其算符、除以 `L` 的位置和返回值；独立 iLTRG 示例、DMRG、实时 TDVP 不因位于同一仓库就自动具有上述契约。`GetResult` 的 `M` 字段只在相应非零场 Chi 分支中赋值，不能假定所有调用都会返回它。

## iLTRG 显式物理接口

`ManyBodySolver/iLTRG/iLTRGForward.m` 使用公共函数
`UtilityFunc/ConvertThermalObservables.m`；常数来自 `QSpinThermalUnits.m`，
与 `GetUnitCon.m` 共用。`result.raw` 同时明确命名自然量和物理量：
`C_natural`、`spin_projection`、`Cm_J_per_mol_K`、`M_muB_per_spin`、
`Chi_emu_per_mol_Oe` 和 `Chi_qspintn_native`。
其中 `Chi_qspintn_native = 4*pi*Chi_emu_per_mol_Oe`，既有契约不变。
`result.data` 的 χ 已是实验 CGS 单位，不能再次按 native → CGS 除以 4π。
该转换目前要求主轴磁场；零场 M/H 未定义。`GetiLTRGRslt`/`iTEBD`
仍返回自然量，需按实际算符和磁场约定转换。

## 逐候选能标与单位边界

多体求解器可共用一个归一化规则：先将 factor、绑定和比例参数展开成以 K 表示的
物理相互作用系数 c，再把所有计算磁场投影到 g 主轴，取
`E0=max(max(abs(c)), max((muB/kB)*norm(g.*B)))`。同一候选的各条曲线共用 E0，
仅相互作用与有效磁场全部为零时用 1 K 并记录原因。此数值尺度不是 CW 实测能标，
也不是全系统谱宽，不额外乘自旋、键数或体系尺寸。

这是 autoQMagen 拟合适配层采用的规则；直接调用旧 GetModel 仍须核对原生能标约定，
直接 iLTRGForward 仍由调用者提供 `energy_scale_K` 和物理 `request.tau`。
相互作用及场能量除以 E0 后，以同一 E0 处理温度和数值步长。ED 自然单位
`beta=E0/T_K`，而 iLTRGForward 的 `raw.beta` 已是 K⁻¹，只取 `T_K=1/raw.beta`，
不能再次乘 E0。公开物理 API 的 `request.tau=tau_dimensionless/E0` K⁻¹；
拟合中应使用本候选解析后的实际 tau，不能直接复制无量纲上限。

例如 iLTRG 拟合默认上限 0.01、E0=2 K 时覆盖至 50 K；需要 100 K 时，autoQMagen
按默认裕量 0.99 取实际 tau=0.00495，物理步长 0.002475 K⁻¹。适配器核对低温
演化步数，求解器仍需 D/tau 收敛验证。XTRG/tanTRG 各自的网格与步长规则不能用
这条 iLTRG 中心差分覆盖公式替代。源数据单位换算、能标归一化、损失归一化应分别记录。

## 常数、温度与场

复现包内结果时使用 `GetUnitCon.m` 的舍入常数：

```text
kB* = 8.617e-2 * 1.602e-22 J/K
muB = 9.274e-24 J/T
NA  = 6.022e23 mol^-1
R   = 8.31446 J/(mol K)
mu0 = 4*pi*1e-7 H/m
h_con = ES*kB*/muB
h = |B|/h_con
T_K = ES/beta
```

`R` 是独立给定的包常数，不用舍入后的 `NA*kB*` 替换。实验磁单位换算依据 [NIST 磁单位表](https://www.nist.gov/document/unitschartpdf)；源码中的量定义和归一化则依据上表列出的函数。

在使用外加场 `B = mu0*H` 的约定下：`H_Oe = 1e4*B_T`，反向 `B_T = 1e-4*H_Oe`，例如 `0.1 T ↔ 1000 Oe`。这里把实验横轴上的 `mu0 H` 转为模型外场；不把材料内部一般关系 `B=mu0(H+M)` 偷换为该外加场约定。记录场方向；横轴若显示 Oe，标作 `H (Oe)`，若显示 T，标明 `B` 或 `mu0 H` 的实际约定。

## 磁性比热 Cm

在自然单位 Hamiltonian 与 `beta` 一致的前提下：

```text
Cm_J_per_mol_K = C_total_natural/L * R
Cm_mJ_per_mol_K = 1000 * Cm_J_per_mol_K
Cm_J_per_mol_K = 0.001 * Cm_mJ_per_mol_K
```

例如 `10 mJ/(mol K) = 0.01 J/(mol K)`。热容中的能量尺度与温度尺度相消，完成上述转换后不要再乘 `ES^2`。若输入本来就是 `Cm_l`，只执行目标显示单位的换算。区分磁性 `Cm` 与含声子背景的总热容 `Cp`；单位换算不能完成背景分离。

## 磁化 M 与 M–H 曲线

沿主轴 α、均匀对角 g 因子的例子，定义

```text
m = g_alpha * <sum_i S_i^alpha>/L       # 数值单位 μB/spin
```

若传入 `ED_chi` 的总磁矩算符已是 `g_alpha*sum_i S_i^alpha`，直接取返回 `M/L`；若只传了总自旋，还需乘 g。**不要重复乘 g，也不要把已经每自旋归一化的输出再次除以 L。**一般 g 张量或不同位点的 g 应由实际磁矩算符处理，不能套用一个任意标量 g。

autoQMagen 默认用 `μB/磁性自旋` 展示磁矩，这是应用层显示约定，显式指定的目标单位优先。明确它与严格 `μB/mol` 的区别。按每摩尔磁性自旋，正反转换为：

| 目标量 | 从 m 转换 | 转回 m |
| --- | --- | --- |
| `M_SI`，A·m²/mol | `NA*muB*m = 5.5848028*m` | `M_SI / 5.5848028` |
| `M_CGS`，emu/mol | `1000*NA*muB*m = 5584.8028*m` | `M_CGS / 5584.8028` |
| 严格 `M_muB_mol`，μB/mol | `NA*m` | `M_muB_mol / NA` |

例如 `S=1/2, gx=1.5`，沿 x 的正向无限强场饱和值为：

```text
m_sat = gx*S = 0.75 μB/spin
      ↔ 4.5165e23 μB/mol
      ↔ 4188.6021 emu/mol
      ↔ 4.1886021 A·m²/mol
```

这是极限值，不保证有限的 8 T 已严格饱和。不能把纵轴 `0.75` 标成 `μB/mol`。每条 M–H 曲线同时记录固定温度、扫场方向、场单位、g 因子和磁矩基准；多个单位列是同一数据的等价表达，不应重复加入损失。

零场磁化可由完整热平均求出，但 `M/H` 在零场不能直接相除。底层 `ED_chi` 在求 M 时不使用参数 `h`，只有其 ratio 输出才除以 `h`；若仅取 M，可用非零 `h` 占位并丢弃 ratio，而 Hamiltonian 必须仍使用真实场。共享 `GetField` 在构造方向时会除以 `norm(B)`，零场需要入口对应的处理。

## 磁化率 χ：原生、CGS 与 SI

`qspintn_native_chi` 是使用方对 `Chi_l` 的标识，不是 MATLAB 字段名。当前代码使用

```text
Chi_con = NA*mu0*muB^2/(ES*kB*) * 1e6
chi_native = (m/h)*Chi_con
           = 1e6*mu0*NA*muB*m/|B|
```

此处 m 已含 g；物理 `Chi_l` 不再额外乘 g²。以下数值换算须保持相同摩尔基准：

| 转换 | 因子 |
| --- | --- |
| 原生 SI cm³/mol → CGS emu/(mol Oe) | 除以 `4*pi` |
| CGS emu/(mol Oe) → 原生 SI cm³/mol | 乘以 `4*pi` |
| 原生 SI cm³/mol → SI m³/mol | 乘以 `1e-6` |
| SI m³/mol → 原生 SI cm³/mol | 乘以 `1e6` |
| CGS emu/(mol Oe) → SI m³/mol | 乘以 `4*pi*1e-6`；反向除以该因子 |

**CGS cm³/mol 与 CGS emu/(mol Oe) 数值相同**。推荐实验图使用后者以明确单位；不能将不带 SI/CGS 说明的 `cm³/mol` 自动当成原生单位。[NIST 表中的摩尔磁化率换算](https://www.nist.gov/document/unitschartpdf)支持上述 CGS/SI 因子。

从同方向的磁化计算有限场 χ 时，非零外加场下：

```text
chi_CGS_emu_per_mol_Oe = M_CGS_emu_per_mol / (1e4*B_T)
chi_SI_m3_per_mol = mu0*M_SI_A_m2_per_mol / B_T
```

`GetResult(..., 'Chi')` 对应有限场 `M/H`，不是自动取 `B→0` 的微分磁化率，也不含退磁修正。对比实验低场 χ 前检查线性区。横向磁矩若不与 Hamiltonian 对易，需在真实场 Hamiltonian 的本征态中计算完整热平均；不能把 `beta*方差(M)` 当作一般量子静态磁化率，也不能仅按 g² 缩放另一方向的曲线。

## 摩尔基准与数据交接

单位来源与归一化来源分开记录。上述物理热力学输出以每摩尔磁性自旋为基准；实验数据可能按每摩尔化学式、每克或整个样品报告。

- 已知每个化学式含 `n` 个参与模型的磁性自旋时，`Y_per_mol_formula = n*Y_per_mol_spin`，反向除以 n。链长 `L` 不是 n；只有明确改报“每摩尔 L 自旋簇”时才按 L 换算。
- 对质量归一化的 Cm、χ、磁矩，先用明确化学式的摩尔质量 `M_formula`（g/mol）：`Y_per_mol_formula = Y_per_g*M_formula`，再按 n 转为每摩尔磁性自旋。反向操作取倒数。没有明确基准时，不从链长或文件名猜测物质的量。
- autoQMagen 在未声明实验摩尔基准时采用 `per_mole_magnetic_spins` 是**应用层默认约定**，不是 QSpinTN 自动识别了实验归一化。使用该默认时记录假设；实验文件或用户明确的基准优先。
- 保留原始列和单位，注明转换因子、常数版本、输出层级、摩尔基准、场方向与 `M/H` 定义。换算后的数据、图、manifest 和报告保持一致；单位换算不修复杂质、偏置或模型不匹配。
