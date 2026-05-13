%% test_quadratic_mixed_manual
% 手动测试二次内部元 + 线性边界元的原始-混合格式
% 不使用边界加密，独立边界网格步长 = 1.25 * h_in

clear; clc; close all;

% ===================== 1. 定义问题 =====================
geom = Rectg(0, 0, 1, 1);                % 几何
u_exact = @(x,y) sin(pi*x) .* sin(pi*y); % 精确解
f = @(x,y) 2*pi^2 * sin(pi*x) .* sin(pi*y); % 右端项
g = @(x,y) 0;                            % 齐次 Dirichlet

% 边界法向导数精确值 ξ = -∂u/∂n
xi_exact_fun = @(x,y) ...
    (abs(x)<1e-12) * ( -pi * sin(pi*y) ) + ...
    (abs(x-1)<1e-12) * ( -pi * sin(pi*y) ) + ...
    (abs(y)<1e-12) * ( -pi * sin(pi*x) ) + ...
    (abs(y-1)<1e-12) * ( -pi * sin(pi*x) );

% ===================== 2. 网格参数 =====================
h_in = 0.1;                % 内部网格初始尺寸
h_bd = 1.25 * h_in;        % 独立边界网格步长（1.25倍）

% ===================== 3. 生成内部网格（线性）并转换为二次 =====================
[p0, ~, t0_lin] = initmesh(geom, 'Hmax', h_in);
p0 = p0';                  % np×2
t0_lin = t0_lin(1:3,:)';   % nt×3

% 转换为二次网格（6节点）
[p, t] = linear2quadMesh(p0, t0_lin);   % p: np×2, t: nt×6

% 提取内部网格的边界边（基于顶点，用于可能的后处理）
T_aux = auxstructure(t(:,1:3));
e_int = T_aux.bdEdge';     % 2×ne_int

% ===================== 4. 生成独立边界网格 =====================
e_b_nodes = boundaryDivide(geom, h_bd);   % 2×nb
nb = size(e_b_nodes, 2);
% 构造闭合边连接（线性边界元，每个节点一个自由度）
e_boundary = [1:nb-1; 2:nb];
e_boundary = [e_boundary, [nb; 1]];       % 2×ne，且 ne = nb

% ===================== 5. 组装刚度矩阵 A =====================
% 注意：StiffnessAssembler2D 期望 p 为 2×np，t 为 6×nt（二次元）
p_for_assembly = p';     % 2×np
t_for_assembly = t';     % 6×nt
A = StiffnessAssembler2D(p_for_assembly, t_for_assembly, @(x,y) 1, 'quadratic');

% ===================== 6. 组装载荷向量 F =====================
F = LoadAssembler2D(p_for_assembly, t_for_assembly, f, 'quadratic');

% ===================== 7. 组装边界耦合矩阵 B 和右端项 G =====================
% BoundaryCoupling2D_quadratic 期望 p 为 np×2，t 为 nt×6
[B, G] = BoundaryCoupling2D_quadratic(p, t, e_b_nodes, e_boundary, g);

% ===================== 8. 求解鞍点系统 =====================
np = size(p, 1);
ne_b = size(e_boundary, 2);
K = [A, B; B', sparse(ne_b, ne_b)];
rhs = [F; G];
sol = K \ rhs;

u = sol(1:np);
xi = -sol(np+1:end);   % 实际法向导数

u = full(u);
xi = full(xi);

% ===================== 9. 误差计算 =====================
% 数值解误差（在节点上）
p_plot = p;   % np×2
U_exact_node = u_exact(p_plot(:,1), p_plot(:,2));
err_u = norm(U_exact_node - u) / sqrt(np);
fprintf('===== 二次内部元 + 线性边界元 =====\n');
fprintf('节点数: %d, 单元数: %d\n', np, size(t,1));
fprintf('数值解相对 L2 误差: %e\n', err_u);

% 边界通量误差（在边界节点上）
xi_exact_node = zeros(nb, 1);
for i = 1:nb
    x = e_b_nodes(1, i);
    y = e_b_nodes(2, i);
    xi_exact_node(i) = xi_exact_fun(x, y);
end
err_xi = norm(xi - xi_exact_node) / sqrt(nb);
fprintf('边界通量相对 L2 误差: %e\n', err_xi);

% ===================== 10. 可视化 =====================
% 数值解云图（注意 trisurf 需要 3×nt 的顶点连接）
figure('Name', '二次元数值解', 'Position', [100, 100, 800, 600]);
trisurf(t(:,1:3), p(:,1), p(:,2), u, 'EdgeColor', 'none');
view(2); axis equal; colorbar; title('u_h (二次元)'); xlabel('x'); ylabel('y');

% 精确解云图
U_exact_plot = u_exact(p(:,1), p(:,2));
figure('Name', '精确解', 'Position', [100, 100, 800, 600]);
trisurf(t(:,1:3), p(:,1), p(:,2), U_exact_plot, 'EdgeColor', 'none');
view(2); axis equal; colorbar; title('u_{exact}'); xlabel('x'); ylabel('y');

% 误差云图
error_plot = U_exact_plot - u;
figure('Name', '误差云图', 'Position', [100, 100, 800, 600]);
trisurf(t(:,1:3), p(:,1), p(:,2), error_plot, 'EdgeColor', 'none');
view(2); axis equal; colorbar; title('误差 u - u_h'); xlabel('x'); ylabel('y');

% 边界通量对比（沿边界节点序号）
s = 1:nb;
figure('Name', '边界通量对比', 'Position', [100, 100, 800, 600]);
plot(s, xi, 'b-o', 'LineWidth', 1.5); hold on;
plot(s, xi_exact_node, 'r--x', 'LineWidth', 1.5);
xlabel('边界节点序号'); ylabel('ξ'); legend('ξ_h', 'ξ_{exact}');
title('边界通量对比（二次内部元+线性边界元）'); grid on;