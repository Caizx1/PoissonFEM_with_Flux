%% test_compare_uniform_vs_concentrated.m
% 对比线性元与二次元在均匀网格和边界加密网格下的收敛性
% 边界加密时内部网格在边界附近单元尺寸约为 h_in^2，独立边界网格步长也取 h_bd = h_in^2

clear; clc; close all;

% ===================== 1. 问题定义 =====================
geom = Rectg(0, 0, 1, 1);                % 几何
u_exact = @(x,y) sin(pi*x) .* sin(pi*y);
grad_u_exact = @(x,y) [pi*cos(pi*x).*sin(pi*y), pi*sin(pi*x).*cos(pi*y)];
f = @(x,y) 2*pi^2 * sin(pi*x) .* sin(pi*y);
g = @(x,y) 0;

xi_exact_fun = @(x,y) ...                % 精确边界法向导数 ξ = -∂u/∂n
    (abs(x)<1e-12) * (-pi * sin(pi*y)) + ...
    (abs(x-1)<1e-12) * (-pi * sin(pi*y)) + ...
    (abs(y)<1e-12) * (-pi * sin(pi*x)) + ...
    (abs(y-1)<1e-12) * (-pi * sin(pi*x));

% ===================== 2. 参数设置 =====================
h_list = [0.2, 0.1, 0.05, 0.04];       % 内部网格初始尺寸
ratio_uniform = 1.25;                   % 均匀网格时边界/内部步长比（线性元和二次元通用）
refine_opts_uniform.use = false;        % 均匀网格
refine_opts_uniform.C = 1.5;            % 占位
refine_opts_uniform.max_iter = 10;

refine_opts_concentrated.use = true;    % 边界加密
refine_opts_concentrated.C = 1.5;       % 加密判据常数
refine_opts_concentrated.max_iter = 10;

% 存储误差
% 线性元
err_u_lin_uniform = []; err_xi_lin_uniform = [];
err_u_lin_concentrated = []; err_xi_lin_concentrated = [];
% 二次元
err_u_quad_uniform = []; err_xi_quad_uniform = [];
err_u_quad_concentrated = []; err_xi_quad_concentrated = [];

% ===================== 3. 循环计算 =====================
for h_in = h_list
    fprintf('\n========== h_in = %g ==========\n', h_in);
    
    % ---- 线性元 ----
    % 均匀网格
    h_bd_uniform = ratio_uniform * h_in;
    [u_lin_uni, xi_lin_uni, p_lin_uni, t_lin_uni, ~, e_b_nodes_uni, e_boundary_uni] = ...
        primal_mixed_solver2D(geom, f, g, h_in, h_bd_uniform, 'linear', refine_opts_uniform);
    p_lin_uni = p_lin_uni';   % np×2
    t_lin_uni = t_lin_uni';   % nt×3
    err_u_lin_uni = H1_error(p_lin_uni, t_lin_uni, u_lin_uni, u_exact, grad_u_exact, 'linear');
    err_xi_lin_uni = L2_boundary_error(e_b_nodes_uni, e_boundary_uni, xi_lin_uni, xi_exact_fun);
    err_u_lin_uniform = [err_u_lin_uniform, err_u_lin_uni];
    err_xi_lin_uniform = [err_xi_lin_uniform, err_xi_lin_uni];
    
    % 边界加密网格
    h_bd_concentrated = 2 * h_in^2;   % 与加密后的边界单元尺寸相当
    [u_lin_con, xi_lin_con, p_lin_con, t_lin_con, ~, e_b_nodes_con, e_boundary_con] = ...
        primal_mixed_solver2D(geom, f, g, h_in, h_bd_concentrated, 'linear', refine_opts_concentrated);
    p_lin_con = p_lin_con';
    t_lin_con = t_lin_con';
    err_u_lin_con = H1_error(p_lin_con, t_lin_con, u_lin_con, u_exact, grad_u_exact, 'linear');
    err_xi_lin_con = L2_boundary_error(e_b_nodes_con, e_boundary_con, xi_lin_con, xi_exact_fun);
    err_u_lin_concentrated = [err_u_lin_concentrated, err_u_lin_con];
    err_xi_lin_concentrated = [err_xi_lin_concentrated, err_xi_lin_con];
    
    fprintf('线性元: 均匀网格 dof=%d, err_u=%.3e, err_xi=%.3e\n', size(p_lin_uni,1), err_u_lin_uni, err_xi_lin_uni);
    fprintf('线性元: 加密网格 dof=%d, err_u=%.3e, err_xi=%.3e\n', size(p_lin_con,1), err_u_lin_con, err_xi_lin_con);
    
    % ---- 二次元 ----
    % 均匀网格（使用改进后的边界耦合，允许任意比值）
    [u_quad_uni, xi_quad_uni, p_quad_uni, t_quad_uni, ~, e_b_nodes_uni, e_boundary_uni] = ...
        primal_mixed_solver2D(geom, f, g, h_in, h_bd_uniform, 'quadratic', refine_opts_uniform);
    p_quad_uni = p_quad_uni';   % np×2
    t_quad_uni = t_quad_uni';   % nt×6
    err_u_quad_uni = H1_error(p_quad_uni, t_quad_uni, u_quad_uni, u_exact, grad_u_exact, 'quadratic');
    err_xi_quad_uni = L2_boundary_error(e_b_nodes_uni, e_boundary_uni, xi_quad_uni, xi_exact_fun);
    err_u_quad_uniform = [err_u_quad_uniform, err_u_quad_uni];
    err_xi_quad_uniform = [err_xi_quad_uniform, err_xi_quad_uni];
    
    % 边界加密网格
    [u_quad_con, xi_quad_con, p_quad_con, t_quad_con, ~, e_b_nodes_con, e_boundary_con] = ...
        primal_mixed_solver2D(geom, f, g, h_in, h_bd_concentrated, 'quadratic', refine_opts_concentrated);
    p_quad_con = p_quad_con';
    t_quad_con = t_quad_con';
    err_u_quad_con = H1_error(p_quad_con, t_quad_con, u_quad_con, u_exact, grad_u_exact, 'quadratic');
    err_xi_quad_con = L2_boundary_error(e_b_nodes_con, e_boundary_con, xi_quad_con, xi_exact_fun);
    err_u_quad_concentrated = [err_u_quad_concentrated, err_u_quad_con];
    err_xi_quad_concentrated = [err_xi_quad_concentrated, err_xi_quad_con];
    
    fprintf('二次元: 均匀网格 dof=%d, err_u=%.3e, err_xi=%.3e\n', size(p_quad_uni,1), err_u_quad_uni, err_xi_quad_uni);
    fprintf('二次元: 加密网格 dof=%d, err_u=%.3e, err_xi=%.3e\n', size(p_quad_con,1), err_u_quad_con, err_xi_quad_con);
end

% ===================== 4. 拟合收敛阶 =====================
ord_u_lin_uni = polyfit(log(h_list), log(err_u_lin_uniform), 1);
ord_xi_lin_uni = polyfit(log(h_list), log(err_xi_lin_uniform), 1);
ord_u_lin_con = polyfit(log(h_list), log(err_u_lin_concentrated), 1);
ord_xi_lin_con = polyfit(log(h_list), log(err_xi_lin_concentrated), 1);
ord_u_quad_uni = polyfit(log(h_list), log(err_u_quad_uniform), 1);
ord_xi_quad_uni = polyfit(log(h_list), log(err_xi_quad_uniform), 1);
ord_u_quad_con = polyfit(log(h_list), log(err_u_quad_concentrated), 1);
ord_xi_quad_con = polyfit(log(h_list), log(err_xi_quad_concentrated), 1);

fprintf('\n===== 收敛阶 (基于 h_in) =====\n');
fprintf('线性元 均匀网格: u H1 = %.2f, ξ L2 = %.2f\n', ord_u_lin_uni(1), ord_xi_lin_uni(1));
fprintf('线性元 加密网格: u H1 = %.2f, ξ L2 = %.2f\n', ord_u_lin_con(1), ord_xi_lin_con(1));
fprintf('二次元 均匀网格: u H1 = %.2f, ξ L2 = %.2f\n', ord_u_quad_uni(1), ord_xi_quad_uni(1));
fprintf('二次元 加密网格: u H1 = %.2f, ξ L2 = %.2f\n', ord_u_quad_con(1), ord_xi_quad_con(1));

% ===================== 5. 绘图 =====================
figure('Name', 'H1误差收敛曲线', 'Position', [100, 100, 800, 600]);
loglog(h_list, err_u_lin_uniform, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
loglog(h_list, err_u_lin_concentrated, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 8);
loglog(h_list, err_u_quad_uniform, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 8);
loglog(h_list, err_u_quad_concentrated, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('内部网格初始尺寸 h'); ylabel('H^1 误差');
title('H^1 误差收敛对比');
legend('线性元 均匀', '线性元 加密', '二次元 均匀', '二次元 加密', 'Location', 'best');
grid on;

figure('Name', '边界L2误差收敛曲线', 'Position', [100, 100, 800, 600]);
loglog(h_list, err_xi_lin_uniform, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
loglog(h_list, err_xi_lin_concentrated, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 8);
loglog(h_list, err_xi_quad_uniform, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 8);
loglog(h_list, err_xi_quad_concentrated, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('内部网格初始尺寸 h'); ylabel('边界 L^2 误差');
title('边界通量 L^2 误差收敛对比');
legend('线性元 均匀', '线性元 加密', '二次元 均匀', '二次元 加密', 'Location', 'best');
grid on;

% ===================== 6. 辅助函数（如果尚未定义） =====================
% 这里假设 H1_error 和 L2_boundary_error 已在路径中，否则需要添加定义
% 若不存在，可复制之前的定义。