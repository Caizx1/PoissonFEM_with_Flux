%% test_ratio_convergence.m
% 使用 H1_error 和 L2_boundary_error 函数，研究不同边界网格步长比下二次内部元+间断线性边界元的收敛阶
% 输出每个 ratio 下 u_H1 误差和 xi_L2 误差的收敛阶，并绘制收敛阶随 ratio 变化的曲线

clear; clc; close all;

% ===================== 1. 问题定义 =====================
examples = examples();

ex = 2;
geom   = examples(ex).geom;
u_exact= examples(ex).u_exact;
grad_u_exact = examples(ex).grad_u_exact;
f      = examples(ex).f;
g      = examples(ex).g;
xi_exact_fun = examples(ex).xi_exact_fun;

% ===================== 2. 参数设置 =====================
h_list = [0.2, 0.1, 0.05, 0.025];       % 内部网格尺寸（粗到细）
ratio_list = [1.25, 1.5,1.75 2]; % 边界/内部步长比

refine_opts.use = false;                % 不使用边界加密

% 预分配存储
conv_u = zeros(length(ratio_list), 1);   % u_H1 收敛阶
conv_xi = zeros(length(ratio_list), 1);  % xi_L2 收敛阶

% ===================== 3. 循环计算 =====================
for r = 1:length(ratio_list)
    ratio = ratio_list(r);
    fprintf('\n========== ratio = %g ==========\n', ratio);
    
    err_u = [];   % 存储不同 h_in 下的 H1 误差
    err_xi = [];  % 存储不同 h_in 下的边界 L2 误差
    h_vals = [];  % 存储实际使用的 h_in（用于拟合）
    
    for h_in = h_list
        h_bd = ratio * h_in;
        
        % 调用求解器（二次元）
        result_quad = primal_mixed_solver2D(geom, f, g, h_in, h_bd, 'quadratic', refine_opts);
        u = result_quad.sol.u;
        xi = result_quad.sol.xi;
        p = result_quad.mesh.p';   % np×2
        t = result_quad.mesh.t';   % nt×6
        
        % 计算 u 的 H1 误差（需要梯度函数）
        err_u_val = H1_error(p, t, u, u_exact, grad_u_exact, 'quadratic');
        err_u = [err_u, err_u_val];
        
        % 计算 xi 的 L2 边界误差（使用独立边界网格信息）
        % 注意：e_boundary 是 2×ne，e_b_nodes 是 2×nb

        xi_nodes = result_quad.boundary.xi_nodes;
        xi_exact_disc = zeros(size(xi));
        for i = 1:length(xi)
            xi_exact_disc(i) = xi_exact_fun(xi_nodes(1,i), xi_nodes(2,i));
        end
        err_xi_val = norm(xi - xi_exact_disc) / sqrt(length(xi));

        err_xi = [err_xi, err_xi_val];
        
        h_vals = [h_vals, h_in];
        dof = size(p,1);
        fprintf('h_in = %g, dof = %d, err_u = %e, err_xi = %e\n', ...
                h_in, dof, err_u_val, err_xi_val);
    end
    
    % 拟合收敛阶（基于 log(h) 和 log(err)）
    p_u = polyfit(log(h_vals), log(err_u), 1);
    p_xi = polyfit(log(h_vals), log(err_xi), 1);
    conv_u(r) = p_u(1);
    conv_xi(r) = p_xi(1);
    fprintf('收敛阶: u_H1 = %.2f, xi_L2 = %.2f\n', conv_u(r), conv_xi(r));
end

% ===================== 4. 输出结果 =====================
% 表格
T = table(ratio_list(:), conv_u(:), conv_xi(:), ...
    'VariableNames', {'Ratio', 'Conv_u_H1', 'Conv_xi_L2'});
disp('========== 不同比值下的收敛阶 ==========');
disp(T);

% 绘制收敛阶随 ratio 变化曲线
figure('Name', '收敛阶 vs h_{bd}/h_{in}', 'Position', [100, 100, 800, 600]);
plot(ratio_list, conv_u, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(ratio_list, conv_xi, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('h_{bd} / h_{in}');
ylabel('收敛阶');
title('二次内部元 + 间断线性边界元：不同边界网格步长比下的收敛阶');
legend('u (H^1 误差)', 'ξ (L^2 边界误差)', 'Location', 'best');
grid on;