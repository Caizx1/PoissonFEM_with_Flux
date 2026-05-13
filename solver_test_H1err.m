%% compare_linear_quadratic_uniform_with_norms
% 对比线性元+常数边界元与二次元+间断线性边界元在均匀网格下的收敛性
% 使用 H1 范数误差和边界 L2 范数误差

clear; clc; close all;

% ===================== 1. 定义问题 =====================
examples = examples();

ex = 6;
geom   = examples(ex).geom;
u_exact= examples(ex).u_exact;
grad_u_exact = examples(ex).grad_u_exact;
f      = examples(ex).f;
g      = examples(ex).g;
xi_exact_fun = examples(ex).xi_exact_fun;

% ===================== 2. 网格参数 =====================
h_list = [0.2, 0.1, 0.05, 0.025,0.0125];   % 内部网格尺寸;
h_bd_ratio = 1.25;                  % 独立边界网格步长 = h_in * ratio
refine_opts.use = false;             % 不使用边界集中加密

% 存储误差
err_u_lin_H1 = []; err_xi_lin_L2 = []; dof_lin = [];
err_u_quad_H1 = []; err_xi_quad_L2 = []; dof_quad = [];

% ===================== 3. 循环计算 =====================
k = 0;
for h_in = h_list
    k = k + 1;
    fprintf('\n========== h_in = %g ==========\n', h_in);
    h_bd_ratio = 1.25;
    h_bd = h_bd_ratio * h_in;

    % ---- 线性元 + 常数边界元 ----
    result_lin = primal_mixed_solver2D(geom, f, g, h_in, h_bd, 'linear', refine_opts);
    u_lin  = result_lin.sol.u;
    xi_lin = result_lin.sol.xi;
    p_lin  = result_lin.mesh.p';   % np×2
    t_lin  = result_lin.mesh.t';   % nt×3
    % H1 误差
    err_u_lin_H1_node = H1_error(p_lin, t_lin, u_lin, u_exact, grad_u_exact, 'linear');
    err_u_lin_H1 = [err_u_lin_H1, err_u_lin_H1_node];
    dof_lin = [dof_lin, size(p_lin,1)];
    % 边界通量 L2 误差（常数边界元：xi_nodes 为边中点）
    xi_nodes_lin = result_lin.boundary.xi_nodes;
    xi_exact = zeros(size(xi_lin));
    for i = 1:length(xi_lin)
        xi_exact(i) = xi_exact_fun(xi_nodes_lin(1,i), xi_nodes_lin(2,i));
    end
    err_xi_lin_L2_node = norm(xi_lin - xi_exact) / sqrt(length(xi_lin));
    err_xi_lin_L2 = [err_xi_lin_L2, err_xi_lin_L2_node];
    fprintf('线性元: dof = %d, H1 err_u = %e, L2 err_xi = %e\n', size(p_lin,1), err_u_lin_H1_node, err_xi_lin_L2_node);

    if k > 1
        order_u_lin = log(err_u_lin_H1(k-1) / err_u_lin_H1(k)) / log(h_list(k-1) / h_list(k));
        order_xi_lin = log(err_xi_lin_L2(k-1) / err_xi_lin_L2(k)) / log(h_list(k-1) / h_list(k));
        fprintf('u收敛阶 = %.2f，ξ收敛阶 = %.2f \n',order_u_lin,order_xi_lin);
    end

    % ---- 二次元 + 间断线性边界元 ----
    h_bd_ratio = 1.25;
    h_bd = h_bd_ratio * h_in;
    result_quad = primal_mixed_solver2D(geom, f, g, h_in, h_bd, 'quadratic', refine_opts);
    u_quad  = result_quad.sol.u;
    xi_quad = result_quad.sol.xi;
    p_quad  = result_quad.mesh.p';   % np×2
    t_quad  = result_quad.mesh.t';   % nt×6
    % H1 误差（二次元）
    err_u_quad_H1_node = H1_error(p_quad, t_quad, u_quad, u_exact, grad_u_exact, 'quadratic');
    err_u_quad_H1 = [err_u_quad_H1, err_u_quad_H1_node];
    dof_quad = [dof_quad, size(p_quad,1)];
    % 边界通量 L2 误差（间断线性边界元：xi_nodes 为高斯点）
    xi_nodes_quad = result_quad.boundary.xi_nodes;
    xi_exact_disc = zeros(size(xi_quad));
    for i = 1:length(xi_quad)
        xi_exact_disc(i) = xi_exact_fun(xi_nodes_quad(1,i), xi_nodes_quad(2,i));
    end
    err_xi_quad_L2_node = norm(xi_quad - xi_exact_disc) / sqrt(length(xi_quad));
    err_xi_quad_L2 = [err_xi_quad_L2, err_xi_quad_L2_node];
    fprintf('二次元: dof = %d, H1 err_u = %e, L2 err_xi = %e\n', size(p_quad,1), err_u_quad_H1_node, err_xi_quad_L2_node);

    if k > 1
        order_u_quad = log(err_u_quad_H1(k-1) / err_u_quad_H1(k)) / log(h_list(k-1) / h_list(k));
        order_xi_quad = log(err_xi_quad_L2(k-1) / err_xi_quad_L2(k)) / log(h_list(k-1) / h_list(k));
        fprintf('u收敛阶 = %.2f，ξ收敛阶 = %.2f \n',order_u_quad,order_xi_quad);
    end
end

% ===================== 4. 绘制收敛曲线 =====================
figure('Name', '收敛曲线对比 (H1 norm)', 'Position', [100, 100, 800, 600]);
loglog(dof_lin, err_u_lin_H1, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
loglog(dof_quad, err_u_quad_H1, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('自由度 (节点数)'); ylabel('相对 H1 误差');
title('H1 范数误差对比');
legend('线性元', '二次元', 'Location', 'southwest');
grid on;

figure('Name', '收敛曲线对比 (L2 boundary)', 'Position', [100, 100, 800, 600]);
loglog(dof_lin, err_xi_lin_L2, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
loglog(dof_quad, err_xi_quad_L2, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('自由度 (节点数)'); ylabel('边界 L2 误差');
title('边界通量 L2 误差对比');
legend('线性元', '二次元', 'Location', 'southwest');
grid on;

% 拟合收敛阶
ord_u_lin = polyfit(log(h_list), log(err_u_lin_H1), 1);
ord_xi_lin = polyfit(log(h_list), log(err_xi_lin_L2), 1);
ord_u_quad = polyfit(log(h_list), log(err_u_quad_H1), 1);
ord_xi_quad = polyfit(log(h_list), log(err_xi_quad_L2), 1);
fprintf('\n线性元收敛阶 (u H1): %.2f\n', ord_u_lin(1));
fprintf('线性元收敛阶 (ξ L2): %.2f\n', ord_xi_lin(1));
fprintf('二次元收敛阶 (u H1): %.2f\n', ord_u_quad(1));
fprintf('二次元收敛阶 (ξ L2): %.2f\n', ord_xi_quad(1));

%% ===================== 辅助函数 =====================
function err = H1_error(p, t, u_h, u_exact, grad_u_exact, order)
    % p: np×2 节点坐标
    % t: nt×3 (linear) 或 nt×6 (quadratic)
    % u_h: 节点解 (np×1)
    % u_exact: 精确解函数句柄
    % grad_u_exact: 精确梯度函数句柄，返回 [∂u/∂x, ∂u/∂y]
    % order: 'linear' 或 'quadratic'

    nt = size(t,1);
    err = 0;

    % 高斯积分点（面积坐标）和权重（三点积分，精确到5次）
    gp = [1/6, 1/6, 2/3, 1/6;
          1/6, 2/3, 1/6, 1/6;
          2/3, 1/6, 1/6, 1/6];
    nGauss = size(gp,1);

    for K = 1:nt
        if strcmp(order, 'linear')
            nodes = t(K,1:3);
            nNodes = 3;
        else
            nodes = t(K,:);
            nNodes = 6;
        end
        xy = p(nodes, :);   % nNodes × 2
        ue = u_h(nodes);    % nNodes × 1

        errK = 0;
        for i = 1:nGauss
            L1 = gp(i,1); L2 = gp(i,2); L3 = gp(i,3); w = gp(i,4);
            if strcmp(order, 'linear')
                % 线性形函数（面积坐标）
                N = [L1, L2, L3];   % 1×3
                % 形函数对 L1, L2 的导数（对 L3 = 1-L1-L2）
                dNdL = [1, 0; 0, 1; -1, -1];  % 3×2
            else
                % 二次形函数
                N = [L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1), 4*L1*L2, 4*L2*L3, 4*L3*L1]; % 1×6
                % 对 L1, L2 的偏导数（6×2）
                dNdL = zeros(6,2);
                dNdL(1,1) = 4*L1 - 1;  dNdL(1,2) = 0;
                dNdL(2,1) = 0;         dNdL(2,2) = 4*L2 - 1;
                dNdL(3,1) = 1 - 4*L3;  dNdL(3,2) = 1 - 4*L3;
                dNdL(4,1) = 4*L2;      dNdL(4,2) = 4*L1;
                dNdL(5,1) = -4*L2;     dNdL(5,2) = 4*L3 - 4*L2;
                dNdL(6,1) = 4*L3 - 4*L1; dNdL(6,2) = -4*L1;
            end
            % Jacobian 矩阵 (2×2)
            J = xy' * dNdL;   % xy' 是 2×nNodes，dNdL 是 nNodes×2，乘积 2×2
            detJ = abs(det(J));
            % 物理坐标
            xq = N * xy(:,1);   % 点积：1×nNodes * nNodes×1 = 标量
            yq = N * xy(:,2);
            % 形函数对物理坐标的导数 (nNodes×2)
            dNdx = dNdL / J;
            % 数值解和梯度
            u_h_q = N * ue;          % 1×nNodes * nNodes×1 = 标量
            grad_u_h_q = ue' * dNdx;  % 1×nNodes * nNodes×2 = 1×2
            % 精确解和梯度
            u_exact_q = u_exact(xq, yq);
            grad_u_exact_q = grad_u_exact(xq, yq);
            % 误差
            e = u_exact_q - u_h_q;
            grad_e = grad_u_exact_q - grad_u_h_q;
            % 累加
            errK = errK + (e^2 + grad_e(1)^2 + grad_e(2)^2) * w * detJ;
        end
        err = err + errK;
    end
    err = sqrt(err);
end
