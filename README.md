# Project 1 in C — 星系团热气体静力学平衡 + 铁的湍流扩散

Computational Astrophysics AA 2026-27，第一次作业的完整 C 实现（Part I + Part II）。
原作业给出的是 Fortran 版本；本目录是逐项重写、结构化、带自检的 C11 版本，
不依赖任何外部库（只需要 `libm`）。

- **Part I**：NFW 暗物质 + Hernquist 中心椭圆星系（BCG）+ 热气体静力学平衡；
  等温、非等温 `T(r)`、以及额外的湍流非热压强。
- **Part II**：ICM 中铁的湍流扩散 PDE，含 SNIa + 恒星风源项；复刻讲义的
  三个实验：只扩散、只加源、源 + 扩散。

---

## 1. 编译

需要 C11 编译器和 `make`：

```sh
gcc --version         # gcc / clang / MinGW-w64 均可
make
```

Windows（MinGW-w64 / w64devkit / MSYS2）同样直接 `make`。如果不想用 make：

```sh
gcc -O2 -std=c11 -Wall -Wextra -o project1 src/main.c src/cluster.c src/diffusion.c -lm
```

## 2. 运行

```sh
./project1 --check     # 先跑自检：15 项数值检查，全部 PASS 才继续
./project1             # 完整计算，结果写入 output/
./project1 --conv      # 网格收敛表（等温模型 vs 解析解）
./project1 --help
```

常用选项：

| 选项 | 含义 | 默认 |
|---|---|---|
| `-n, --jmax N` | 网格边界数 | 5000 |
| `-R, --rmax KPC` | 外半径 | 2800 kpc |
| `--rvir KPC` | 计算重子分数的半径 | 2800 kpc |
| `--r500 KPC` | `T(r)` 中的 R500 | 1400 kpc |
| `--fb X` | 目标重子分数 | 0.16 |
| `--rho0 X` | 固定中心气体密度（跳过求解） | 自动求解 |
| `--vt KM/S` | 湍流速度 | 500 |
| `--turb-iters N` | 湍流不动点最大迭代数 | 100 |
| `--snu X` | SNIa 率（SNu） | 0.15 |
| `--kappa X` | 扩散系数 [cm²/s] | 1.32e29 |
| `--kappa-vary` | 用 `D(r)`（中心降低 60%） | 关 |
| `--alpha-evolve` | 用含时 `α(t)` 的源项 | 关 |
| `--no-fe` | 只算 Part I | 关 |
| `-o, --out DIR` | 输出目录 | `output` |

例：

```sh
./project1 --snu 0.72                 # 讲义 p.46 的“铁够了”模型
./project1 --kappa 1.06e30 --snu 2.5  # 讲义 p.47 的“铁够了但分布太集中”
./project1 --rho0 1.6e-25             # 用讲义自己的 rho0 复核
./project1 -n 20000 --no-fe           # 高分辨率 Part I
```

## 3. 输出

`output/` 下都是两列以上的纯文本表，第一行是 `#` 注释，任何绘图工具都能读：

| 文件 | 列 |
|---|---|
| `grid.dat` | `r_left[kpc]  r_center[kpc]` |
| `mass.dat` | `r[kpc]  M_DM_num[Msun]  M_DM_analytic[Msun]` |
| `gas.dat` | `r_c[kpc]  rho_iso  rho_analytic  rho_NFW  [g/cm³]` |
| `gas_BCG.dat` | `r_c[kpc]  rho_iso_BCG  rho_NFW` |
| `gas_BCG_non-iso.dat` | `r_c[kpc]  rho_T(r)  rho_NFW` |
| `temperature.dat` | `r_c[kpc]  T[K]  T/Tmg` |
| `gas_BCG_non-iso-turb.dat` | `r_c[kpc]  rho_turb  rho_noniso  P_nt  P_nt/P_th` |
| `fe_observed.dat` | `r_c[kpc]  Z_obs/Zsun  (Z_obs-Z_out)/Zsun` |
| `fe_diffusion.dat` | `r_c[kpc]  Z(t=0,1,2,5 Gyr)` |
| `fe_source.dat` | `r_c[kpc]  Z(t=1,5 Gyr)` |
| `fe_both.dat` | `r_c[kpc]  Z(t=1,5 Gyr)` |
| `fe_masses.dat` | `t[Gyr]  experiment  M_Fe(<100kpc)[Msun]  M_Fe(tot)[Msun]` |
| `summary.txt` | 一次运行的全部关键数字 |

## 4. 代码结构

```
src/cluster.h     常量、Grid/Cluster 结构、Part I 接口
src/cluster.c     网格、NFW 质量、BCG、T(r)、HSE 积分、rho0 二分、残差
src/diffusion.h   Part II 接口
src/diffusion.c   Fe 的有限体积扩散 + 源项 + 时间步 + 快照
src/main.c        CLI、流程、诊断、表格输出、--check / --conv
tools/plot.mjs    可选的 Node 绘图脚本（读 output/ → figures/）
figures/          本次运行生成的 9 张图
```

关键设计：

- **网格**：`n` 个边界 `r_j`，密度、温度、`Z_Fe` 在壳中心 `r_{j+1/2}`，
  质量、引力在边界；`ΔV_j = 4π/3 (r_{j+1}³ − r_j³)`。
- **第一层暗物质质量**用 NFW 解析质量 `M(r_1)` 代替 `ρ(r_c)ΔV`，
  消除 NFW 中心发散带来的 33% 偏差；其余层仍然数值累加，再与解析质量比较。
- **HSE 积分**用压力形式：`P_{j+1} = P_j − Δr G M_{j+1}/r_{j+1}² (ρ_j+ρ_{j+1})/2`，
  对 `ρ` 隐式、对 `P` 二阶；等价于 `exp(−A/C)` 的 Padé 近似，天然保正。
  非等温时 `P = ρ kT/(μmp)` 直接带 `dT/dr`，不需要单独写温度梯度项。
- **`ρ0`** 不是拍脑袋给的：对 `fb(Rvir)=0.16` 做几何二分（~50 次），
  每次重新积分整条剖面；等温模型还按解析解把 `ρ0` 从 `r=0` 换算到第一层中心。
- **湍流**：`P_tot = P_th + P_nt`，`P_nt = ⅓ρv_t²`。因为 `v_t` 常数时
  方程对 `ρ` 是线性的，固定点迭代收敛很快；程序打印 1–5 次迭代与收敛解的差，
  3 次已经到 `1.6e-4`（讲义写“3 iterations”）。
- **Fe 扩散**：`ρ_Fe` 在壳中心，通量在边界；零通量边界条件 → 总铁质量严格守恒
  （实测漂移 `5e-15`）。时间步 `Δt = 0.4 Δr²/(2D)`，每步重新计算。
- **Fe 源项**：`S = ρ_*(r)[α_*(t) Z_*/1.4 + α_SN(t) M_Fe,Ia/M_Ia]`；
  `ρ_*` 用壳平均 `(M_*(r_{j+1})−M_*(r_j))/ΔV_j`，使体积分严格等于 `M_*`，
  避免 Hernquist 中心 `1/r` 发散造成的几个百分点误差。

## 5. 结果速览

完整输出见 `summary.txt`。`--check` 的关键检查（jmax=5000）：

```
DM mass vs analytic                 7.4e-07
isothermal vs analytic              8.9e-04   (r>10 kpc 最大相对误差)
fb(iso/BCG/noniso/turb)=0.16        4 个都精确到 1e-10
HSE residual, T(r)                  1.2e-03
HSE residual, turbulence            1.1e-03
turbulence: 3 iterations enough     1.6e-04
Fe source total vs analytic         1.8e-16
Fe mass conserved (diffusion)       5.3e-15
```

Part II 与讲义参考值（p.42–46）的对比：

| 实验 | 量 | 本代码 | 讲义 |
|---|---|---|---|
| 只加源，1 Gyr | `M_Fe(<100kpc)` | 1.517e7 | 1.52e7 |
| 只加源，5 Gyr | `M_Fe(<100kpc)` | 7.582e7 | 7.58e7 |
| 源+扩散，1 Gyr | `M_Fe(<100kpc)` | 1.504e7 | 1.50e7 |
| 源+扩散，5 Gyr | `M_Fe(<100kpc)` | 6.423e7 | 6.40e7 |
| 源+扩散，SNu=0.72，5 Gyr | `M_Fe(<100kpc)` / total | 2.805e8 / 3.634e8 | 2.80e8 / 3.63e8 |

（单位都是 `M☉`。只加源的“总量”对 `M_*` 的积分是精确的，所以和解析式
`5.3e-22 M_* t` 相对误差 1.8e-16。）

## 6. 图

`figures/` 里已经有本次运行生成的 9 张 PNG（质量、等温密度、四个密度模型、
`T(r)`、湍流、Fe 三个实验、Fe 质量预算）。重新生成：

```sh
cd tools
npm install            # 只装 @napi-rs/canvas
node plot.mjs
```

如果不想用 Node，`output/*.dat` 可以直接用 gnuplot / Python matplotlib / 任何
工具画；本目录没有把绘图绑进 C 程序。

## 7. 与 Fortran 原版的差异

1. **没有复制粘贴**：4 段 HSE 积分、4 段质量积分合并成 `integrate_model()`；
   4 个模型只差一个枚举。
2. **`ρ0` 是解出来的**，不是事后整体缩放；输出里给出每个模型的 `ρ0` 和 `fb`。
3. **`ρ0` 定义在 `r=0`**（等温模型用解析解换算到第一层中心），
   而不是把中心密度直接放在 `r=0.3 kpc`。
4. **第一层 NFW 质量用解析值**，`mass.dat` 的数值/解析差从 33% 降到 7e-7。
5. **湍流迭代到收敛**，并打印 1–5 次迭代的收敛表；3 次与收敛解差 1.6e-4。
6. **Fe 源项按壳精确积分**，总铁质量与 `5.3e-22 M_* t` 精确到机器精度。
7. 所有输出都有 `--check` 和 `--conv`，不再靠“看图觉得对”。

## 8. 已知局限

- 网格均匀；讲义也允许中心加密，但均匀网格已经是一阶收敛（见 `--conv`）。
- 湍流 `v_t` 常数；`P_nt` 随半径变化的情形只留了 `--kappa-vary` 的扩散系数接口，
  没有做 `v_t(r)`。
- 源项的 `α(t)` 时间演化有实现（`alphaEvolve`），但默认关闭以复刻讲义
  “先当作常数”的参考值；打开它需要改 `diffusion.c` 里的 `fe_source_coef` 调用。
- Part II 的 Fe 只跟气体走，没有考虑其他星系（讲义 p.49 的 addendum）。

## 9. 参考

- 课程讲义 `cas2026-27_proj1.pdf`（Part I p.8–22，Part II p.23–52）。
- Rebusco et al. 2005, 2006（Fe 扩散的观测背景）。
- 原 Fortran 代码 `project_1.f90` 和 `Project1-review.md`（同工作区）。
