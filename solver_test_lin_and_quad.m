%% compare_linear_quadratic_uniform
% 对比线性元+常数边界元与二次元+线性边界元在均匀网格下的收敛性
% 精确解: u = sin(pi*x)*sin(pi*y), 齐次 Dirichlet 边界

clear; clc; close all;

% ===================== 1. 定义问题 =====================
geom = Rectg(0, 0, 1, 1);                % 几何
u_exact = @(x,y) sin(pi*x) .* sin(pi*y); % 精确解
f = @(x,y) 2*pi^2 * sin(pi*x) .* sin(pi*y); % 右端项
g = @(x,y) 0;                            % 齐次 Dirichlet

% 边界法向导数精确值 ξ = -∂u/∂n
xi_exact_fun = @(x,y) ...
    (abs(x)<1e-12) * (-pi * sin(pi*y)) + ...
    (abs(x-1)<1e-12) * (-pi * sin(pi*y)) + ...
    (abs(y)<1e-12) * (-pi * sin(pi*x)) + ...
    (abs(y-1)<1e-12) * (-pi * sin(pi*x));

% ===================== 2. 网格参数 =====================
h_list = [0.2, 0.1, 0.05, 0.025];   % 内部网格尺寸
%h_list = [0.2 0.16 0.1 0.08 0.05];
h_bd_ratio = 1.25;                  % 独立边界网格步长 = h_in * ratio
refine_opts.use = false;             % 不使用边界集中加密

% 存储误差和自由度
err_u_lin = []; err_xi_lin = []; dof_lin = [];
err_u_quad = []; err_xi_quad = []; dof_quad = [];

% ===================== 3. 循环计算 =====================
for h_in = h_list
    fprintf('\n========== h_in = %g ==========\n', h_in);
    h_bd = h_bd_ratio * h_in;

    % ---- 线性元 + 常数边界元 ----
    h_bd_ratio = 1.25;
    h_bd = h_bd_ratio * h_in;
    [u_lin, xi_lin, p_lin, t_lin, e_lin, e_b_nodes, e_bd] = ...
        primal_mixed_solver2D(geom, f, g, h_in, h_bd, 'linear', refine_opts);
    np_lin = size(p_lin, 2);
    nt_lin = size(t_lin, 2);
    % 节点误差
    p_lin_plot = p_lin';                % np×2
    U_exact_lin = u_exact(p_lin_plot(:,1), p_lin_plot(:,2));
    err_u_lin_node = norm(U_exact_lin - u_lin) / sqrt(np_lin);
    err_u_lin = [err_u_lin, err_u_lin_node];
    dof_lin = [dof_lin, np_lin];
    % 边界通量误差（常数边界元：每条边中点）
    midpoints = (e_b_nodes(:, e_bd(1,:)) + e_b_nodes(:, e_bd(2,:))) / 2;
    xi_exact = zeros(size(xi_lin));
    for i = 1:length(xi_lin)
        m = midpoints(:,i);
        xi_exact(i) = xi_exact_fun(m(1), m(2));
    end
    err_xi_lin_node = norm(xi_lin - xi_exact) / sqrt(length(xi_lin));
    err_xi_lin = [err_xi_lin, err_xi_lin_node];
    fprintf('线性元: dof = %d, err_u = %e, err_xi = %e\n', np_lin, err_u_lin_node, err_xi_lin_node);

    % ---- 二次元 + 线性边界元 ----
    h_bd_ratio = 1;
    h_bd = h_bd_ratio * h_in;
    [u_quad, xi_quad, p_quad, t_quad, e_quad, e_b_nodes, e_bd] = ...
        primal_mixed_solver2D(geom, f, g, h_in, h_bd, 'quadratic', refine_opts);
    np_quad = size(p_quad, 2);
    nt_quad = size(t_quad, 2);
    % 节点误差
    p_quad_plot = p_quad';              % np×2
    U_exact_quad = u_exact(p_quad_plot(:,1), p_quad_plot(:,2));
    err_u_quad_node = norm(U_exact_quad - u_quad) / sqrt(np_quad);
    err_u_quad = [err_u_quad, err_u_quad_node];
    dof_quad = [dof_quad, np_quad];
    % 边界通量误差（线性边界元：每个节点）
    xi_exact_node = zeros(size(xi_quad));
    for i = 1:length(xi_quad)
        x = e_b_nodes(1, i);
        y = e_b_nodes(2, i);
        xi_exact_node(i) = xi_exact_fun(x, y);
    end
    err_xi_quad_node = norm(xi_quad - xi_exact_node) / sqrt(length(xi_quad));
    err_xi_quad = [err_xi_quad, err_xi_quad_node];
    fprintf('二次元: dof = %d, err_u = %e, err_xi = %e\n', np_quad, err_u_quad_node, err_xi_quad_node);
end

% ===================== 4. 绘制收敛曲线 =====================
figure('Name', '收敛曲线对比', 'Position', [100, 100, 800, 600]);
loglog(dof_lin, err_u_lin, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
loglog(dof_lin, err_xi_lin, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 8);
loglog(dof_quad, err_u_quad, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 8);
loglog(dof_quad, err_xi_quad, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('自由度 (节点数)'); ylabel('相对 L2 误差');
title('收敛性对比 (均匀网格)');
legend('u 线性元', 'ξ 线性元', 'u 二次元', 'ξ 二次元', 'Location', 'southwest');
grid on;

% 拟合收敛阶
ord_u_lin = polyfit(log(h_list), log(err_u_lin), 1);
ord_xi_lin = polyfit(log(h_list), log(err_xi_lin), 1);
ord_u_quad = polyfit(log(h_list), log(err_u_quad), 1);
ord_xi_quad = polyfit(log(h_list), log(err_xi_quad), 1);
fprintf('\n线性元收敛阶 (u): %.2f\n', ord_u_lin(1));
fprintf('线性元收敛阶 (ξ): %.2f\n', ord_xi_lin(1));
fprintf('二次元收敛阶 (u): %.2f\n', ord_u_quad(1));
fprintf('二次元收敛阶 (ξ): %.2f\n', ord_xi_quad(1));

% ===================== 5. 可视化 h_in = 0.1 时的结果 =====================
h_vis = 0.1;
h_bd_ratio = 1.25;
h_bd_vis = h_bd_ratio * h_vis;

% 5.1 线性元结果
[u_lin_vis, xi_lin_vis, p_lin_vis, t_lin_vis, ~, e_b_nodes_lin, e_bd_lin] = ...
    primal_mixed_solver2D(geom, f, g, h_vis, h_bd_vis, 'linear', refine_opts);
p_lin_vis = p_lin_vis';   % np×2
t_lin_vis = t_lin_vis';   % nt×3
U_exact_lin_vis = u_exact(p_lin_vis(:,1), p_lin_vis(:,2));
err_lin_vis = U_exact_lin_vis - u_lin_vis;

% 线性元误差云图
figure('Name', '线性元误差云图 (h=0.1)', 'Position', [100, 100, 800, 600]);
trisurf(t_lin_vis, p_lin_vis(:,1), p_lin_vis(:,2), err_lin_vis, 'EdgeColor', 'none');
view(2); axis equal; colorbar; title('误差 u - u_h (线性元, h=0.1)'); xlabel('x'); ylabel('y');

% 线性元边界通量对比（常数边界元，每条边一个值）
midpoints_lin = (e_b_nodes_lin(:, e_bd_lin(1,:)) + e_b_nodes_lin(:, e_bd_lin(2,:))) / 2;
xi_exact_lin = zeros(size(xi_lin_vis));
for i = 1:length(xi_lin_vis)
    m = midpoints_lin(:,i);
    xi_exact_lin(i) = xi_exact_fun(m(1), m(2));
end
figure('Name', '线性元边界通量对比 (h=0.1)', 'Position', [100, 100, 800, 600]);
s_lin = 1:length(xi_lin_vis);
plot(s_lin, xi_lin_vis, 'b-o', 'LineWidth', 1.5); hold on;
plot(s_lin, xi_exact_lin, 'r--x', 'LineWidth', 1.5);
xlabel('边界边序号'); ylabel('ξ = -∂u/∂n');
legend('ξ_h', 'ξ_{exact}'); title('线性元边界通量对比 (h=0.1)'); grid on;

% 5.2 二次元结果
h_bd_ratio = 1;
h_bd = h_bd_ratio * h_in;
[u_quad_vis, xi_quad_vis, p_quad_vis, t_quad_vis, ~, e_b_nodes_quad, e_bd_quad] = ...
    primal_mixed_solver2D(geom, f, g, h_vis, h_bd_vis, 'quadratic', refine_opts);
p_quad_vis = p_quad_vis';   % np×2
t_quad_vis = t_quad_vis';   % nt×3（仅顶点）
U_exact_quad_vis = u_exact(p_quad_vis(:,1), p_quad_vis(:,2));
err_quad_vis = U_exact_quad_vis - u_quad_vis;

% 二次元误差云图（使用顶点连接）
figure('Name', '二次元误差云图 (h=0.1)', 'Position', [100, 100, 800, 600]);
trisurf(t_quad_vis(:,1:3), p_quad_vis(:,1), p_quad_vis(:,2), err_quad_vis, 'EdgeColor', 'none');
view(2); axis equal; colorbar; title('误差 u - u_h (二次元, h=0.1)'); xlabel('x'); ylabel('y');

% 二次元边界通量对比（线性边界元，每个节点一个值）
xi_exact_quad_node = zeros(size(xi_quad_vis));
for i = 1:length(xi_quad_vis)
    x = e_b_nodes_quad(1, i);
    y = e_b_nodes_quad(2, i);
    xi_exact_quad_node(i) = xi_exact_fun(x, y);
end
figure('Name', '二次元边界通量对比 (h=0.1)', 'Position', [100, 100, 800, 600]);
s_quad = 1:length(xi_quad_vis);
plot(s_quad, xi_quad_vis, 'b-o', 'LineWidth', 1.5); hold on;
plot(s_quad, xi_exact_quad_node, 'r--x', 'LineWidth', 1.5);
xlabel('边界节点序号'); ylabel('ξ = -∂u/∂n');
legend('ξ_h', 'ξ_{exact}'); title('二次元边界通量对比 (h=0.1)'); grid on;