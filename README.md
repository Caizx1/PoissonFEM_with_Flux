# Poisson 方程有限元求解器集合

本项目实现了多种求解 Poisson 方程 & -\Delta u = f & 的有限元方法，支持 **Dirichlet 边界条件**，并可计算边界法向导数 $ \xi = \partial u/\partial n$。包含三种主要求解器：

- **原始-混合格式**（primal-mixed）：同时求解 $u$ 和边界法向导数 $\xi$，支持独立边界网格（非匹配网格），可处理不连续通量。
- **标准有限元 + 直接投影法**：先求解 $u$，再在单元边界上直接计算 $\xi$（通过梯度点乘法向量）。
- **标准有限元 + 变分法（残差法）**：求解 $u$ 后，利用残差和边界质量矩阵计算 $\xi$。

## 主要特性

- 支持 **线性元** 和 **二次元**（六节点三角形）。
- 支持 **均匀网格** 和 **边界集中加密网格**。
- 原始-混合格式支持 **独立边界网格**（边界与内部网格可非匹配）。
- 边界通量可选用 **常数元**、**连续线性元** 或 **间断线性元**（高斯点自由度）。
- 提供多个测试算例（光滑解、角点奇异性、非齐次边界等）。

## 目录结构
`plain
PoissonFEM_with_Flux/
├── solvers/ # 主求解器
│ ├── primal_mixed_solver2D.m # 原始-混合格式
│ ├── fem_normal_derivative.m # 标准有限元 + 直接法向导数
│ └── fem_variational_derivative.m # 标准有限元 + 变分法向导数
│
├── mesh/ # 网格生成与处理
│ ├── linear2quadMesh.m
│ ├── trimeshRefine.m
│ ├── boundary_concentrated_mesh.m
│ ├── ...
│
├── fem_assembly/ # 有限元矩阵组装
│ ├── StiffnessAssembler2D.m
│ ├── LoadAssembler2D.m
│ ├── HatGradients.m
│ └── boundaryMassAssembler2D.m
│
├── boundary_coupling/ # 原始-混合格式的边界耦合（独立网格）
│ ├── BoundaryCoupling2D_independent.m
│ ├── BoundaryCoupling2D_discontinuous_linear.m
│ ├── segment_overlap.m
│ └── BoundaryCoupling2D_quadratic.m # 连续线性边界元
│
├── projection/ # 法向导数投影（后处理）
│ ├── project_normal_derivative_p1.m
│ ├── project_normal_derivative_p2.m
│ └── outer_normal.m
│
├── error_estimation/ # 误差计算
│ ├── H1_error.m
│ ├── L2_boundary_error.m
│ └── L2_error_discontinuous_linear.m
│
├── examples/ # 算例库
│ ├── examples.m
│ └── Rectg.m
│
├── utils/ # 通用辅助函数
│ ├── init_solver_params.m
│ ├── isInTriangle.m
│ ├── triLocate.m
│ ├── auxstructure.m # 来自 iFEM
│ └── bisect.m # 来自 iFEM
│
├── tests/ # 测试脚本
│ ├── test_ratio_convergence.m
│ ├── test_fem_variational_derivative.m
│ ├── test_compare_uniform_vs_concentrated.m
│ └── ... # 更多测试
│
└── README.md
└── setpath.m
└── rm_path.m
`

## 环境依赖

- MATLAB R2019b 或更高版本（推荐）。
- **iFEM 有限元包**（提供 `auxstructure`, `bisect` 等函数，已内置在utils中）。  
  下载地址：https://github.com/lyc102/ifem  
  安装方法：将 iFEM 文件夹加入 MATLAB 路径，或运行其 `setpath.m` 脚本。

## 快速开始

### 1. 添加项目路径
```matlab
% 在项目根目录下运行 setpath 脚本，自动添加所有子文件夹到 MATLAB 路径
setpath
```

### 2. 运行简单示例（齐次边界，光滑解）
```matlab
geom = Rectg(0,0,1,1);
u_exact = @(x,y) sin(pi*x).*sin(pi*y);
f = @(x,y) 2*pi^2 * sin(pi*x).*sin(pi*y);
g = @(x,y) 0;   % 齐次 Dirichlet

% 使用原始-混合格式求解（线性元 + 常数边界元）
result = primal_mixed_solver2D(geom, f, g, 0.1, 0.125, 'linear');
u = result.sol.u;
xi = result.sol.xi;   % 边界法向导数（常数元）

% 绘制云图
p = result.mesh.p'; t = result.mesh.t';
trisurf(t(:,1:3), p(:,1), p(:,2), u, 'EdgeColor', 'none');
```

### 3. 使用标准有限元 + 变分法向导数（二次元）
```matlab
[sol, mesh] = fem_variational_derivative(geom, f, g, 0.1, 'quadratic');
u = sol.u;
lambda = sol.lambda;   % 边界节点上的法向导数
```
## 参考
- PFEFFERER J, WINKLER M. Finite element error estimates for normal derivatives on boundary concentrated meshes[J/OL]. SIAM Journal on Numerical Analysis,2019, 57(5): 2043-2073. DOI: 10.1137/18M1181341.

- GATICA G N. Springerbriefs in mathematics: A simple introduction to the mixed finite element method: Theory and applications[M/OL]. Cham: Springer International Publishing, 2014. DOI: 10.1007/978-3-319-03695-3.

- L. Chen. iFEM: an integrated finite element method package in MATLAB. Technical Report, University of California at Irvine, 2009.