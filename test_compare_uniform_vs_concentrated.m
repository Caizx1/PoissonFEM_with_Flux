%% test_compare_uniform_vs_concentrated.m
% 对比线性元（常数边界元）与二次元（间断线性边界元）
% 在均匀网格和边界加密网格下的收敛性
% 边界加密时内部网格在边界附近单元尺寸约为 h_in^2，独立边界网格步长也取 h_bd = h_in^2

clear; clc; close all;

% ===================== 1. 问题定义 =====================
examples = examples();
ex = 1;   % 可选择其他算例
geom   = examples(ex).geom;
u_exact= examples(ex).u_exact;
grad_u_exact = examples(ex).grad_u_exact;
f      = examples(ex).f;
g      = examples(ex).g;
xi_exact_fun = examples(ex).xi_exact_fun;

% ===================== 2. 参数设置 =====================
h_list = [0.2, 0.1, 0.05, 0.04];       % 内部网格初始尺寸
ratio_uniform = 1.25;                   % 均匀网格时边界/内部步长比
refine_opts_uniform.use = false;        % 均匀网格
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
    
    % ---- 线性元（常数边界元） ----
    h_bd_uniform = ratio_uniform * h_in;
    % 均匀网格
    res_lin_uni = primal_mixed_solver2D(geom, f, g, h_in, h_bd_uniform, 'linear', refine_opts_uniform);
    err_u_lin_uni = H1_error(res_lin_uni.mesh.p', res_lin_uni.mesh.t', ...
                             res_lin_uni.sol.u, u_exact, grad_u_exact, 'linear');
    err_xi_lin_uni = L2_error_constant_boundary(res_lin_uni.boundary.e_b_nodes, ...
                                                res_lin_uni.boundary.e_boundary, ...
                                                res_lin_uni.sol.xi, xi_exact_fun);
    err_u_lin_uniform = [err_u_lin_uniform, err_u_lin_uni];
    err_xi_lin_uniform = [err_xi_lin_uniform, err_xi_lin_uni];
    
    % 加密网格
    h_bd_concentrated = 2 * h_in^2;
    res_lin_con = primal_mixed_solver2D(geom, f, g, h_in, h_bd_concentrated, 'linear', refine_opts_concentrated);
    err_u_lin_con = H1_error(res_lin_con.mesh.p', res_lin_con.mesh.t', ...
                             res_lin_con.sol.u, u_exact, grad_u_exact, 'linear');
    err_xi_lin_con = L2_error_constant_boundary(res_lin_con.boundary.e_b_nodes, ...
                                                res_lin_con.boundary.e_boundary, ...
                                                res_lin_con.sol.xi, xi_exact_fun);
    err_u_lin_concentrated = [err_u_lin_concentrated, err_u_lin_con];
    err_xi_lin_concentrated = [err_xi_lin_concentrated, err_xi_lin_con];
    
    fprintf('线性元: 均匀网格 dof=%d, err_u=%.3e, err_xi=%.3e\n', ...
            size(res_lin_uni.mesh.p,2), err_u_lin_uni, err_xi_lin_uni);
    fprintf('线性元: 加密网格 dof=%d, err_u=%.3e, err_xi=%.3e\n', ...
            size(res_lin_con.mesh.p,2), err_u_lin_con, err_xi_lin_con);
    
    % ---- 二次元（间断线性边界元） ----
    % 均匀网格
    res_quad_uni = primal_mixed_solver2D(geom, f, g, h_in, h_bd_uniform, 'quadratic', refine_opts_uniform);
    err_u_quad_uni = H1_error(res_quad_uni.mesh.p', res_quad_uni.mesh.t', ...
                              res_quad_uni.sol.u, u_exact, grad_u_exact, 'quadratic');
    err_xi_quad_uni = L2_error_discontinuous_linear(res_quad_uni.boundary.e_b_nodes, ...
                                                    res_quad_uni.boundary.e_boundary, ...
                                                    res_quad_uni.sol.xi, xi_exact_fun);
    err_u_quad_uniform = [err_u_quad_uniform, err_u_quad_uni];
    err_xi_quad_uniform = [err_xi_quad_uniform, err_xi_quad_uni];
    
    % 加密网格
    res_quad_con = primal_mixed_solver2D(geom, f, g, h_in, h_bd_concentrated, 'quadratic', refine_opts_concentrated);
    err_u_quad_con = H1_error(res_quad_con.mesh.p', res_quad_con.mesh.t', ...
                              res_quad_con.sol.u, u_exact, grad_u_exact, 'quadratic');
    err_xi_quad_con = L2_error_discontinuous_linear(res_quad_con.boundary.e_b_nodes, ...
                                                    res_quad_con.boundary.e_boundary, ...
                                                    res_quad_con.sol.xi, xi_exact_fun);
    err_u_quad_concentrated = [err_u_quad_concentrated, err_u_quad_con];
    err_xi_quad_concentrated = [err_xi_quad_concentrated, err_xi_quad_con];
    
    fprintf('二次元: 均匀网格 dof=%d, err_u=%.3e, err_xi=%.3e\n', ...
            size(res_quad_uni.mesh.p,2), err_u_quad_uni, err_xi_quad_uni);
    fprintf('二次元: 加密网格 dof=%d, err_u=%.3e, err_xi=%.3e\n', ...
            size(res_quad_con.mesh.p,2), err_u_quad_con, err_xi_quad_con);
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
title('边界通量 L^2 误差对比');
legend('线性元 均匀', '线性元 加密', '二次元 均匀', '二次元 加密', 'Location', 'best');
grid on;

% ===================== 辅助误差函数 =====================
function err = L2_error_constant_boundary(e_b_nodes, e_boundary, xi, xi_exact_fun)
    % e_b_nodes: 2×nb 边界节点坐标
    % e_boundary: 2×ne 边界边连接
    % xi: 长度 ne 的向量，每个边上常数
    % xi_exact_fun: 精确法向导数函数句柄
    ne = size(e_boundary, 2);
    % 高斯积分点（两点，精确到三次）
    gauss_pts = [-1/sqrt(3), 1/sqrt(3)];
    gauss_w = [1, 1];
    err2 = 0;
    for e = 1:ne
        idx1 = e_boundary(1, e);
        idx2 = e_boundary(2, e);
        p1 = e_b_nodes(:, idx1);
        p2 = e_b_nodes(:, idx2);
        len = norm(p2 - p1);
        xi_h = xi(e);
        % 在边上积分 (ξ_h - ξ_exact)^2
        for q = 1:2
            t = (gauss_pts(q) + 1) / 2;   % 映射到 [0,1]
            x = (1-t)*p1(1) + t*p2(1);
            y = (1-t)*p1(2) + t*p2(2);
            xi_ex = xi_exact_fun(x, y);
            err2 = err2 + (xi_h - xi_ex)^2 * (len/2) * gauss_w(q);
        end
    end
    err = sqrt(err2);
end

function err = L2_error_discontinuous_linear(e_b_nodes, e_boundary, xi, xi_exact_fun)
    % e_b_nodes: 2×nb 边界节点坐标
    % e_boundary: 2×ne 边界边连接
    % xi: 长度 2*ne 的向量，每个单元两个自由度（高斯点）
    % xi_exact_fun: 精确法向导数函数句柄
    ne = size(e_boundary, 2);
    % 高斯点参数 t ∈ [0,1]
    t_gauss = [0.5 - 0.5/sqrt(3), 0.5 + 0.5/sqrt(3)];
    % 参考区间 [-1,1] 上的高斯积分（用于单元内积分，精确到三次）
    tau_gauss = [-1/sqrt(3), 1/sqrt(3)];
    w_gauss = [1, 1];
    err2 = 0;
    for e = 1:ne
        idx1 = e_boundary(1, e);
        idx2 = e_boundary(2, e);
        p1 = e_b_nodes(:, idx1);
        p2 = e_b_nodes(:, idx2);
        len = norm(p2 - p1);
        xi1 = xi(2*(e-1)+1);
        xi2 = xi(2*(e-1)+2);
        % 基函数（间断线性）
        psi1 = @(t) (t - t_gauss(2)) / (t_gauss(1) - t_gauss(2));
        psi2 = @(t) (t - t_gauss(1)) / (t_gauss(2) - t_gauss(1));
        % 在边上积分 (ξ_h - ξ_exact)^2
        % 使用两点高斯积分（在 [0,1] 上，映射自 [-1,1]）
        for q = 1:2
            tau = tau_gauss(q);
            t = 0.5 + 0.5*tau;   % 映射到 [0,1]
            x = (1-t)*p1(1) + t*p2(1);
            y = (1-t)*p1(2) + t*p2(2);
            xi_h = psi1(t)*xi1 + psi2(t)*xi2;
            xi_ex = xi_exact_fun(x, y);
            err2 = err2 + (xi_h - xi_ex)^2 * (len/2) * w_gauss(q);
        end
    end
    err = sqrt(err2);
end