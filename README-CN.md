# Mini Synergetics and Dissipative Structures（迷你协同论与耗散结构）

一套**从零实现、零外部依赖的 C 语言**协同学与耗散结构理论实现库——研究远离热力学平衡系统中自组织、序参量涌现与图灵斑图形成。每个模块将 Haken（协同学）、Prigogine（耗散结构）、Turing（形态发生）和 Landau（相变理论）的奠基性工作转化为可运行、自包含的 C 代码。

## 子模块

| 子模块 | 主题 | 核心课程 |
|--------|------|----------|
| [mini-chemical-oscillations-bz-reaction](mini-chemical-oscillations-bz-reaction/) | Belousov-Zhabotinsky 反应，Field-Koros-Noyes (FKN) 机理，Oregonator 模型，质量作用动力学，Hopf 分岔，极限环稳定性，双稳态分析，螺旋波与靶波，反应-扩散耦合，BZ 系统中的 Turing 不稳定性 | MIT 5.60, MIT 18.385J |
| [mini-entropy-production-minimization](mini-entropy-production-minimization/) | 熵产生率 σ，Prigogine 最小熵产生原理 (MEPP)，Onsager 倒易关系，涨落定理（Crooks, Jarzynski），非线性不可逆热力学，Glansdorff-Prigogine 稳定性判据，稳态分析，耗散系统热力学 | MIT 10.40, Stanford ChemE 340 |
| [mini-haken-slaving-principle](mini-haken-slaving-principle/) | Haken 役使原理 (Versklavungsprinzip)，快变量绝热消去，时间尺度层级分解，谱隙条件，序参量识别（谱方法/统计/物理方法），非线性协同学动力学，临界点附近的稳定性与分岔 | Santa Fe Institute CSSS, MIT 18.385J |
| [mini-nonequilibrium-phase-transitions](mini-nonequilibrium-phase-transitions/) | 作为非平衡相变的分岔理论，Landau-Ginzburg 平均场理论，临界现象（幂律奇异性，普适类），序参量与关联函数，随机动力学（主方程，Fokker-Planck, Langevin），远离平衡态的对称性破缺 | MIT 8.334, Cornell PHYS 7654 |
| [mini-order-parameters-emergence](mini-order-parameters-emergence/) | 序参量动力学 (dξ/dt = αξ − βξ³ + 噪声)，涌现量化指标（超额熵，ε-机复杂度），集体动力学（Kuramoto 同步，Vicsek 蜂拥），Haken 役使原理，Landau-Ginzburg 自由能，非平衡相变，ENSO 应用 | Santa Fe Institute CSSS, MIT 8.334 |
| [mini-prigogine-dissipative-structures](mini-prigogine-dissipative-structures/) | Prigogine-Lefever Brusselator 模型，反应-扩散 PDE 系统，线性/非线性稳定性分析，分岔检测（鞍结、叉形、Hopf、跨临界），耗散结构中的熵产生，Onsager 倒易关系，Glansdorff-Prigogine 稳定性判据，扩展不可逆热力学 | MIT 10.40, MIT 5.60 |
| [mini-self-organization-far-equilibrium](mini-self-organization-far-equilibrium/) | 远离平衡态的自组织机制，分岔分析（规范型，余维数），非线性动力学积分，时空斑图分类（条纹、斑点、六边形、螺旋波、行波），非平衡热力学（熵、自由能、化学势），Prigogine-Haken-Nicolis 综合理论 | Santa Fe Institute CSSS, MIT 18.385J |
| [mini-turing-pattern-formation](mini-turing-pattern-formation/) | Turing 不稳定性（活化-抑制机理），7 种经典反应-扩散模型（Gierer-Meinhardt, Gray-Scott, Brusselator, Schnakenberg, FitzHugh-Nagumo, Lengyel-Epstein, Thomas），色散关系 λ(k)，Turing 空间构建，PDE 数值求解器（显式/隐式/RK4），斑图分类（斑点、条纹、迷宫） | MIT 18.354J, MIT 18.385J |

## 设计理念

- **零外部依赖**——纯 C (C99/C11)，仅依赖 `libc` 和 `libm`
- **模块自包含**——每个子目录拥有独立的 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码的映射**——每个结构体和函数均对标 Haken、Prigogine、Turing 或 Landau 的奠基性著作
- **从零实现**——不使用任何现有库；每个概念均从热力学、协同学和反应-扩散基本原理构建

## 构建

每个模块独立运行。进入模块目录后执行：

```bash
cd mini-chemical-oscillations-bz-reaction
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-synergetics-dissipative/
├── mini-chemical-oscillations-bz-reaction/  # BZ 反应，Oregonator，FKN 机理，螺旋波
├── mini-entropy-production-minimization/    # 最小熵产生原理，Onsager 关系，涨落定理
├── mini-haken-slaving-principle/            # 役使原理，绝热消去，序参量
├── mini-nonequilibrium-phase-transitions/   # Landau 理论，临界现象，随机动力学
├── mini-order-parameters-emergence/         # 序参量动力学，涌现量化，集体动力学
├── mini-prigogine-dissipative-structures/   # Brusselator，反应-扩散，分岔，热力学
├── mini-self-organization-far-equilibrium/  # 自组织机制，斑图分类
└── mini-turing-pattern-formation/           # Turing 不稳定性，7 种 RD 模型，PDE 求解器
```

## 许可证

MIT
