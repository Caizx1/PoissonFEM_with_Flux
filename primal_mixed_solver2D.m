function [u,xi,p,t,e,e_b_nodes,e_boundary] = ...
    primal_mixed_solver2D(geom,f,g,h_in,h_bd,order,refine_opts)
% primal_mixed_solver2D，使用原始-混合格式求解Poisson方程
% -Δu = f，in Ω，u = g on ∂Ω
% 返回数值解 u 和边界法向导数 ∂u/∂n
% 输入：
%   geom：几何矩阵
%   f：右端项函数句柄f(x,y)
%   g：Dirichelet边界条件函数句柄g(x,y)
%   h_in：内部网格尺寸
%   h_bd：边界独立网格尺寸
%   order：内部基函数阶数，'linear'（默认）或'quadratic'
%   refine_opts：结构体，控制边界加密
%       .use：是否使用边界加密（默认false不使用）
%       .C：判据常数（默认1.5）
%       .max_iter：边界加密最大迭代次数（默认10）
% 输出：
%   u：数值解（列向量）
%   xi：边界法向导数的近似（列向量）
%   p：内部网格节点坐标（2*np）
%   t：内部网格三角形单元（3*nt）
%   e：内部网格边界边矩阵（2*ne，每一列为端点索引）
%   e_b_nodes：边界独立网格节点坐标（2*ne_b）
%   e_boundary：边界独立网格边信息（2*ne_b，每一列为端点索引）

% 设置默认参数
if nargin < 6
    order = 'linear';
end

if nargin < 7
    refine_opts = struct();
end
if ~isfield(refine_opts, 'use')
    refine_opts.use = false;
end
if ~isfield(refine_opts, 'C')
    refine_opts.C = 1.5;
end
if ~isfield(refine_opts, 'max_iter')
    refine_opts.max_iter = 10;
end

% 内部网格
[p0,~,t0] = initmesh(geom,"Hmax",h_in);
p0 = p0';
t0 = t0(1:3,:)';

% 边界加密
if refine_opts.use
    [p,e,t] = boundary_concentrated_mesh...
        (p0,t0,h_in,geom,refine_opts.C,refine_opts.max_iter);
else
    p = p0;
    t = t0;
    T = auxstructure(t);
    e = T.bdEdge';
end

if strcmp(order,'quadratic')
    [p,t] = linear2quadMesh(p,t);
    T = auxstructure(t(:,1:3));
    e = T.bdEdge';
end

p = p';
t = t';

% 边界独立网格
e_b_nodes = boundaryDivide(geom,h_bd);
ne_node= size(e_b_nodes,2);
e_boundary = [1:ne_node;2:ne_node,1];

% 矩阵组装
if strcmp(order,'linear')
    A = StiffnessAssembler2D(p, t, @(x,y) 1);
    F = LoadAssembler2D(p, t, f);
    [B, G] = BoundaryCoupling2D_independent(p, e, e_b_nodes, g);
else
    A = StiffnessAssembler2D(p,t,@(x,y)1,'quadratic');
    F = LoadAssembler2D(p,t,f,'quadratic');
    [B,G] = BoundaryCoupling2D_discontinuous_linear(p',t',e_b_nodes,e_boundary,g);
end

% 线性系统求解
np = size(p,2);
nb_dof = size(B,2);
K = [A,B;B',sparse(nb_dof,nb_dof)];
rhs = [F;G];
sol = K \ rhs;

u = sol(1:np);
xi = -sol(np + 1:end);

u = full(u);
xi = full(xi);