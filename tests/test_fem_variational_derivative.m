%% test_fem_variational_derivative.m
% 测试变分法（残差+边界质量矩阵）计算法向导数的收敛性
% 对比线性元与二次元

clear; clc; close all;

% ===================== 1. 选择算例 =====================
examples = examples();
ex = 6;   % 选择算例编号
geom = examples(ex).geom;
u_exact = examples(ex).u_exact;
grad_u_exact = examples(ex).grad_u_exact;
f = examples(ex).f;
g = examples(ex).g;
xi_exact_fun = examples(ex).xi_exact_fun;

% ===================== 2. 参数设置 =====================
h_list = [0.2, 0.1, 0.05, 0.025];   % 网格尺寸
refine_opts.use = false;             % 不使用边界加密

% 存储误差
err_u_lin_H1 = []; err_lambda_lin_L2 = []; h_actual_lin = [];
err_u_quad_H1 = []; err_lambda_quad_L2 = []; h_actual_quad = [];

% ===================== 3. 循环计算 =====================
for h = h_list
    fprintf('\n========== h = %g ==========\n', h);
    
    % ---- 线性元 ----
    [sol_lin, mesh_lin] = fem_variational_derivative(geom, f, g, h, 'linear', refine_opts);
    u_lin = sol_lin.u;
    lambda_lin = sol_lin.lambda;          % 边界节点上的法向导数
    p_lin = mesh_lin.p;
    t_lin = mesh_lin.t;
    e_lin = mesh_lin.e;                   % 边界边（顶点对）
    lambda_loc_lin = mesh_lin.lambda_location; % 边界节点坐标 (n_bd × 2)
    
    % 实际平均边界边长（用于横坐标）
    edge_lengths = zeros(size(e_lin,1),1);
    for i = 1:size(e_lin,1)
        v1 = e_lin(i,1); v2 = e_lin(i,2);
        edge_lengths(i) = norm(p_lin(v1,:) - p_lin(v2,:));
    end
    h_actual_lin = [h_actual_lin, mean(edge_lengths)];
    
    % H1 误差
    err_u_lin = H1_error(p_lin, t_lin, u_lin, u_exact, grad_u_exact, 'linear');
    err_u_lin_H1 = [err_u_lin_H1, err_u_lin];
    
    % 边界 L2 误差（加权质量矩阵）
    err_lambda_lin = L2_error_boundary_weighted(p_lin, e_lin, lambda_loc_lin, lambda_lin, xi_exact_fun, 'linear');
    err_lambda_lin_L2 = [err_lambda_lin_L2, err_lambda_lin];
    
    fprintf('线性元: dof=%d, H1 err_u=%.3e, L2 err_lambda=%.3e\n', size(p_lin,1), err_u_lin, err_lambda_lin);
    
    % ---- 二次元 ----
    [sol_quad, mesh_quad] = fem_variational_derivative(geom, f, g, h, 'quadratic', refine_opts);
    u_quad = sol_quad.u;
    lambda_quad = sol_quad.lambda;
    p_quad = mesh_quad.p;
    t_quad = mesh_quad.t;
    e_quad = mesh_quad.e;
    lambda_loc_quad = mesh_quad.lambda_location; % n_bd × 2
    % 平均边界边长
    edge_lengths = zeros(size(e_quad,1),1);
    for i = 1:size(e_quad,1)
        v1 = e_quad(i,1); v2 = e_quad(i,2);
        edge_lengths(i) = norm(p_quad(v1,:) - p_quad(v2,:));
    end
    h_actual_quad = [h_actual_quad, mean(edge_lengths)];
    
    % H1 误差
    err_u_quad = H1_error(p_quad, t_quad, u_quad, u_exact, grad_u_exact, 'quadratic');
    err_u_quad_H1 = [err_u_quad_H1, err_u_quad];
    
    % 边界 L2 误差（加权质量矩阵）
    err_lambda_quad = L2_error_boundary_weighted(p_quad, e_quad, lambda_loc_quad, lambda_quad, xi_exact_fun, 'quadratic');
    err_lambda_quad_L2 = [err_lambda_quad_L2, err_lambda_quad];
    
    fprintf('二次元: dof=%d, H1 err_u=%.3e, L2 err_lambda=%.3e\n', size(p_quad,1), err_u_quad, err_lambda_quad);
end

% ===================== 4. 收敛阶拟合 =====================
ord_u_lin = polyfit(log(h_actual_lin), log(err_u_lin_H1), 1);
ord_lambda_lin = polyfit(log(h_actual_lin), log(err_lambda_lin_L2), 1);
ord_u_quad = polyfit(log(h_actual_quad), log(err_u_quad_H1), 1);
ord_lambda_quad = polyfit(log(h_actual_quad), log(err_lambda_quad_L2), 1);

fprintf('\n===== 收敛阶 (基于平均边界边长) =====\n');
fprintf('线性元: u H1 = %.2f, λ L2 = %.2f\n', ord_u_lin(1), ord_lambda_lin(1));
fprintf('二次元: u H1 = %.2f, λ L2 = %.2f\n', ord_u_quad(1), ord_lambda_quad(1));

% ===================== 5. 绘图 =====================
figure('Name', 'H1误差收敛', 'Position', [100, 100, 800, 600]);
loglog(h_actual_lin, err_u_lin_H1, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
loglog(h_actual_quad, err_u_quad_H1, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('平均边界边长 h'); ylabel('H^1 误差');
title('H^1 误差收敛对比 (变分法法向导数)');
legend('线性元', '二次元', 'Location', 'southwest');
grid on;

figure('Name', '边界L2误差收敛', 'Position', [100, 100, 800, 600]);
loglog(h_actual_lin, err_lambda_lin_L2, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
loglog(h_actual_quad, err_lambda_quad_L2, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('平均边界边长 h'); ylabel('边界 L^2 误差');
title('法向导数 L^2 误差对比 (变分法)');
legend('线性元', '二次元', 'Location', 'southwest');
grid on;

% ===================== 辅助误差函数 =====================
function err = L2_error_boundary_weighted(p, e, lambda_loc, lambda, xi_exact_fun, order)
% 使用边界质量矩阵计算加权 L2 误差
% p: np×2, e: ne×2 (边界边顶点对)
% lambda_loc: n_bd×2，边界节点坐标
% lambda: n_bd×1，边界节点上的法向导数值
% xi_exact_fun: 精确法向导数函数句柄
% order: 'linear' 或 'quadratic'
    if strcmp(order, 'linear')
        % 线性元：边界节点为顶点
        bdNodes = unique(e(:));
        n_bd = length(bdNodes);
        node2idx = zeros(size(p,1),1);
        node2idx(bdNodes) = 1:n_bd;
        M_b = sparse(n_bd, n_bd);
        for i = 1:size(e,1)
            n1 = e(i,1); n2 = e(i,2);
            idx1 = node2idx(n1); idx2 = node2idx(n2);
            p1 = p(n1,:); p2 = p(n2,:);
            L = norm(p2-p1);
            Me = [2,1;1,2] * L / 6;
            M_b([idx1,idx2],[idx1,idx2]) = M_b([idx1,idx2],[idx1,idx2]) + Me;
        end
        exact = zeros(n_bd,1);
        for i = 1:n_bd
            exact(i) = xi_exact_fun(lambda_loc(i,1), lambda_loc(i,2));
        end
        err_vec = lambda - exact;
        err = sqrt(err_vec' * M_b * err_vec);
    else
        % 二次元：边界节点包括顶点和中点
        n_bd = size(lambda_loc,1);
        ne = size(e,1);
        % 构建每条边的三个节点在 lambda_loc 中的索引
        edge_nodes = zeros(ne,3); % [v1, mid, v2]
        tol = 1e-10;
        for i = 1:ne
            v1 = e(i,1); v2 = e(i,2);
            p1 = p(v1,:); p2 = p(v2,:);
            mid_coord = (p1 + p2) / 2;
            % 查找索引
            idx1 = find(all(abs(lambda_loc - p1) < tol, 2), 1);
            idx2 = find(all(abs(lambda_loc - p2) < tol, 2), 1);
            idx3 = find(all(abs(lambda_loc - mid_coord) < tol, 2), 1);
            if isempty(idx1) || isempty(idx2) || isempty(idx3)
                error('无法找到边界节点坐标匹配');
            end
            edge_nodes(i,:) = [idx1, idx3, idx2];
        end
        % 组装边界质量矩阵
        M_b = sparse(n_bd, n_bd);
        for i = 1:ne
            idx1 = edge_nodes(i,1); idx2 = edge_nodes(i,3); idx3 = edge_nodes(i,2);
            p1 = p(e(i,1),:); p2 = p(e(i,2),:);
            L = norm(p2-p1);
            Me = [4, -1, 2; -1, 4, 2; 2, 2, 16] * L / 30;
            M_b([idx1,idx2,idx3],[idx1,idx2,idx3]) = M_b([idx1,idx2,idx3],[idx1,idx2,idx3]) + Me;
        end
        % 精确值
        exact = zeros(n_bd,1);
        for i = 1:n_bd
            exact(i) = xi_exact_fun(lambda_loc(i,1), lambda_loc(i,2));
        end
        err_vec = lambda - exact;
        err = sqrt(err_vec' * M_b * err_vec);
    end
end

function M = boundary_mass_matrix_linear(p, e, bdNodes)
% 线性元边界质量矩阵（用于误差计算时的辅助）
% 此处返回完整矩阵（bdNodes 应为所有边界节点索引）
    ne = size(e,1);
    M = sparse(length(bdNodes), length(bdNodes));
    node2idx = zeros(size(p,1),1);
    node2idx(bdNodes) = 1:length(bdNodes);
    for i = 1:ne
        n1 = e(i,1); n2 = e(i,2);
        idx1 = node2idx(n1); idx2 = node2idx(n2);
        p1 = p(n1,:); p2 = p(n2,:);
        L = norm(p2-p1);
        Me = [2,1;1,2] * L / 6;
        M([idx1,idx2],[idx1,idx2]) = M([idx1,idx2],[idx1,idx2]) + Me;
    end
end