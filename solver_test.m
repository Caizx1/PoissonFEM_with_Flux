%% test_primal_mixed_solver2D_nonhomogeneous
% 测试非齐次 Dirichlet 边界条件的 Poisson 方程
% 精确解: u(x,y) = x^2 + y^2
% 边界条件: g = x^2 + y^2
% 右端项: f = -4

clear; clc; close all;

% ===================== 定义问题和精确解 =====================
geom = Rectg(0, 0, 1, 1);   % 几何

u_exact = @(x,y) sin(pi * x) .* sin(pi * y);
f = @(x,y) 2 * pi ^ 2 * sin(pi * x) .* sin(pi * y);
g = @(x,y) 0;
xi_exact_fun = @(x,y) ...
    (abs(x)<1e-12) * (-pi) * sin(pi * y) + ...
    (abs(x-1)<1e-12) * (-pi) * sin(pi * y) + ...
    (abs(y)<1e-12) * (-pi) * sin(pi * x) + ...
    (abs(y-1)<1e-12) * (-pi) * sin(pi * x);

% ===================== 网格参数设置 =====================
h_list = [0.2, 0.1, 0.05, 0.025];   % 内部网格尺寸
h_bd_ratio = 1.25;   % 独立边界网格步长 = h_in * ratio，保证稳定性
refine_opts.C = 1.5;
refine_opts.max_iter = 10;

% 存储误差
err_u_uniform = []; err_xi_uniform = [];
err_u_refined = []; err_xi_refined = [];
dof_uniform = []; dof_refined = [];

for h_in = h_list
    fprintf('\n========== h_in = %g ==========\n', h_in);
    h_bd = h_bd_ratio * h_in;   % 满足 h_bd > h_in
    
    % ---- 均匀网格（不加密） ----
    refine_opts.use = false;
    [u_unif, xi_unif, p_unif, t_unif, e_unif, e_b_nodes, e_bd] = ...
        primal_mixed_solver2D(geom, f, g, h_in, h_bd,'linear',refine_opts);
    
    % 计算节点误差（离散 L2 范数）
    p_np2 = p_unif';          % np×2
    u_exact_node = u_exact(p_np2(:,1), p_np2(:,2));
    err_u_unif_node = norm(u_exact_node - u_unif) / sqrt(size(p_unif,2));
    err_u_uniform = [err_u_uniform, err_u_unif_node];
    dof_uniform = [dof_uniform, size(p_unif,2)];
    
    % 边界通量误差（在独立边界边中点处）
    midpoints = (e_b_nodes(:, e_bd(1,:)) + e_b_nodes(:, e_bd(2,:))) / 2;
    xi_exact_vec = zeros(size(xi_unif));

    for i = 1:length(xi_unif)
        m = midpoints(:,i);
        xi_exact_vec(i) = xi_exact_fun(m(1), m(2));
    end

    % figure('Name', '边界通量对比', 'Position', [100, 100, 800, 600]);
    % s = 1:length(xi_unif);
    % plot(s, xi_unif, 'b-o', 'LineWidth', 1.5); hold on;
    % plot(s, xi_exact_vec, 'r--x', 'LineWidth', 1.5);
    % xlabel('边界边序号'); ylabel('ξ = -∂u/∂n');
    % legend('ξ_h', 'ξ_{exact}'); title('边界通量对比');
    % grid on;

    err_xi_unif_node = norm(xi_unif - xi_exact_vec) / sqrt(length(xi_unif));
    err_xi_uniform = [err_xi_uniform, err_xi_unif_node];
    
    fprintf('均匀网格: dof = %d, err_u = %e, err_xi = %e\n', ...
        size(p_unif,2), err_u_unif_node, err_xi_unif_node);
    
    % ---- 边界集中加密网格 ----
    refine_opts.use = true;
    % h_bd = 2 * h_in ^2;
    [u_ref, xi_ref, p_ref, t_ref, e_ref, e_b_nodes, e_bd] = ...
        primal_mixed_solver2D(geom, f, g, h_in, h_bd,'linear', refine_opts);
    
    p_np2_ref = p_ref';
    u_exact_node_ref = u_exact(p_np2_ref(:,1), p_np2_ref(:,2));
    err_u_ref_node = norm(u_exact_node_ref - u_ref) / sqrt(size(p_ref,2));
    err_u_refined = [err_u_refined, err_u_ref_node];
    dof_refined = [dof_refined, size(p_ref,2)];
    
    midpoints_ref = (e_b_nodes(:, e_bd(1,:)) + e_b_nodes(:, e_bd(2,:))) / 2;
    xi_exact_vec_ref = zeros(size(xi_ref));
    for i = 1:length(xi_ref)
        m = midpoints_ref(:,i);
        xi_exact_vec_ref(i) = xi_exact_fun(m(1), m(2));
    end
    err_xi_ref_node = norm(xi_ref - xi_exact_vec_ref) / sqrt(length(xi_ref));
    err_xi_refined = [err_xi_refined, err_xi_ref_node];
    
    fprintf('加密网格: dof = %d, err_u = %e, err_xi = %e\n', ...
        size(p_ref,2), err_u_ref_node, err_xi_ref_node);
end

% ===================== 绘制收敛曲线 =====================
figure('Name', '收敛曲线', 'Position', [100, 100, 800, 600]);
loglog(dof_uniform, err_u_uniform, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
loglog(dof_uniform, err_xi_uniform, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 8);
loglog(dof_refined, err_u_refined, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 8);
loglog(dof_refined, err_xi_refined, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('自由度 (节点数)'); ylabel('相对 L2 误差');
title('收敛性对比');
legend('u 均匀网格', 'ξ 均匀网格', 'u 边界集中网格', 'ξ 边界集中网格', ...
    'Location', 'southwest');
grid on;

% 添加参考线（收敛阶）
ref_h = h_list(end);
ref_dof = dof_uniform(end);
ref_err = err_u_uniform(end);
% 拟合收敛阶（简单示例）
ord_u_uniform = polyfit(log(h_list), log(err_u_uniform), 1);
ord_xi_uniform = polyfit(log(h_list),log(err_xi_uniform),1);
ord_u_refined = polyfit(log(h_list), log(err_u_refined), 1);
ord_xi_refined = polyfit(log(h_list),log(err_xi_refined),1);
fprintf('\n均匀网格收敛阶 (u): %.2f\n', ord_u_uniform(1));
fprintf('均匀网格收敛阶(xi)：%.2f\n',ord_xi_uniform(1));
fprintf('边界集中网格收敛阶 (u): %.2f\n', ord_u_refined(1));
fprintf('边界集中网格收敛阶(xi)：%.2f\n',ord_xi_refined(1));

% ===================== 可视化一个中等网格的结果 =====================
h_vis = 0.1;
h_bd_vis = h_bd_ratio * h_vis;
refine_opts.use = true;   % 使用加密网格展示
[u_vis, xi_vis, p_vis, t_vis, e_vis, e_b_nodes, e_bd] = ...
    primal_mixed_solver2D(geom, f, g, h_vis, h_bd_vis,'linear', refine_opts);

p_plot = p_vis'; t_plot = t_vis';
U_exact_plot = u_exact(p_plot(:,1), p_plot(:,2));

% 数值解云图
figure('Name', '数值解 u_h', 'Position', [100, 100, 800, 600]);
trisurf(t_plot, p_plot(:,1), p_plot(:,2), u_vis, 'EdgeColor', 'none');
view(3); axis equal; colorbar; title('数值解 u_h'); xlabel('x'); ylabel('y');

% 精确解云图
figure('Name', '精确解 u_{exact}', 'Position', [100, 100, 800, 600]);
trisurf(t_plot, p_plot(:,1), p_plot(:,2), U_exact_plot, 'EdgeColor', 'none');
view(3); axis equal; colorbar; title('精确解 u_{exact}'); xlabel('x'); ylabel('y');

% 误差云图
error_plot = U_exact_plot - u_vis;
figure('Name', '误差 u - u_h', 'Position', [100, 100, 800, 600]);
trisurf(t_plot, p_plot(:,1), p_plot(:,2), error_plot, 'EdgeColor', 'none');
view(2); axis equal; colorbar; title('误差 u - u_h'); xlabel('x'); ylabel('y');

% 边界通量对比
midpoints_vis = (e_b_nodes(:, e_bd(1,:)) + e_b_nodes(:, e_bd(2,:))) / 2;
xi_exact_vis = zeros(size(xi_vis));
for i = 1:length(xi_vis)
    m = midpoints_vis(:,i);
    xi_exact_vis(i) = xi_exact_fun(m(1), m(2));
end
figure('Name', '边界通量对比', 'Position', [100, 100, 800, 600]);
s = 1:length(xi_vis);
plot(s, xi_vis, 'b-o', 'LineWidth', 1.5); hold on;
plot(s, xi_exact_vis, 'r--x', 'LineWidth', 1.5);
xlabel('边界边序号'); ylabel('ξ = -∂u/∂n');
legend('ξ_h', 'ξ_{exact}'); title('边界通量对比');
grid on;

% 绘制独立边界网格和内部网格边界
% figure('Name', '网格边界贴合情况', 'Position', [100, 100, 800, 600]);
% triplot(t_plot, p_plot(:,1), p_plot(:,2), 'k-', 'LineWidth', 0.5); hold on;
% plot(e_b_nodes(1, e_bd), e_b_nodes(2, e_bd), 'r-', 'LineWidth', 2);
% axis equal; title('内部网格边界（黑）与独立边界网格（红）');