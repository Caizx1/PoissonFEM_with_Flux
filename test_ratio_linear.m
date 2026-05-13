%% test_ratio_linear.m
% 测试线性内部元 + 常数边界元在不同边界网格步长比下的收敛阶
% 要求 ratio > 1（边界网格比内部网格粗）
% 使用 H1_error 和自定义的常数边界 L2 误差函数

clear; clc; close all;

% ===================== 1. 问题定义 =====================
geom = Rectg(0, 0, 1, 1);                % 几何
u_exact = @(x,y) sin(pi*x) .* sin(pi*y);
grad_u_exact = @(x,y) [pi*cos(pi*x).*sin(pi*y), pi*sin(pi*x).*cos(pi*y)];
f = @(x,y) 2*pi^2 * sin(pi*x) .* sin(pi*y);
g = @(x,y) 0;                            % 齐次 Dirichlet

xi_exact_fun = @(x,y) ...                % 精确边界法向导数 ξ = -∂u/∂n
    (abs(x)<1e-12) * (-pi * sin(pi*y)) + ...
    (abs(x-1)<1e-12) * (-pi * sin(pi*y)) + ...
    (abs(y)<1e-12) * (-pi * sin(pi*x)) + ...
    (abs(y-1)<1e-12) * (-pi * sin(pi*x));

% ===================== 2. 参数设置 =====================
h_list = [0.2, 0.1, 0.05, 0.025];       % 内部网格尺寸
ratio_list = [1.25,1.5,1.6,1.75, 2, 2.5, 3];          % 边界/内部步长比，严格大于 1

refine_opts.use = false;                % 不使用边界加密

% 预分配存储
conv_u = zeros(length(ratio_list), 1);   % u_H1 收敛阶
conv_xi = zeros(length(ratio_list), 1);  % xi_L2 收敛阶

% ===================== 3. 循环计算 =====================
for r = 1:length(ratio_list)
    ratio = ratio_list(r);
    fprintf('\n========== ratio = %g ==========\n', ratio);
    
    err_u = [];
    err_xi = [];
    h_vals = [];
    
    for h_in = h_list
        h_bd = ratio * h_in;
        
        % 调用求解器（线性元）
        [u, xi, p, t, e, e_b_nodes, e_boundary] = ...
            primal_mixed_solver2D(geom, f, g, h_in, h_bd, 'linear', refine_opts);
        
        p = p';               % 2×np → np×2
        t = t';               % 4×nt → nt×3 (线性元)
        
        % 计算 u 的 H1 误差
        err_u_val = H1_error(p, t, u, u_exact, grad_u_exact, 'linear');
        err_u = [err_u, err_u_val];
        
        % 计算 xi 的 L2 边界误差（常数边界元：每条边一个值）
        % 注意：对于常数边界元，xi 长度为 ne (边数)，对应每条边
        % 使用边中点采样
        ne = size(e_boundary, 2);
        xi_exact_mid = zeros(ne, 1);
        for i = 1:ne
            idx1 = e_boundary(1,i);
            idx2 = e_boundary(2,i);
            p1 = e_b_nodes(:,idx1);
            p2 = e_b_nodes(:,idx2);
            mid = (p1 + p2) / 2;
            xi_exact_mid(i) = xi_exact_fun(mid(1), mid(2));
        end
        err_xi_val = norm(xi - xi_exact_mid) / sqrt(ne);
        err_xi = [err_xi, err_xi_val];
        
        h_vals = [h_vals, h_in];
        dof = size(p,1);
        fprintf('h_in = %g, dof = %d, err_u = %e, err_xi = %e\n', ...
                h_in, dof, err_u_val, err_xi_val);
    end
    
    % 拟合收敛阶（基于 h）
    p_u = polyfit(log(h_vals), log(err_u), 1);
    p_xi = polyfit(log(h_vals), log(err_xi), 1);
    conv_u(r) = p_u(1);
    conv_xi(r) = p_xi(1);
    fprintf('收敛阶: u_H1 = %.2f, xi_L2 = %.2f\n', conv_u(r), conv_xi(r));
end

% ===================== 4. 输出结果 =====================
T = table(ratio_list(:), conv_u(:), conv_xi(:), ...
    'VariableNames', {'Ratio', 'Conv_u_H1', 'Conv_xi_L2'});
disp('========== 不同比值下的收敛阶 ==========');
disp(T);

% 绘制收敛阶随 ratio 变化曲线
figure('Name', '线性元：收敛阶 vs h_{bd}/h_{in}', 'Position', [100, 100, 800, 600]);
plot(ratio_list, conv_u, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(ratio_list, conv_xi, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('h_{bd} / h_{in}');
ylabel('收敛阶');
title('线性内部元 + 常数边界元：不同边界网格步长比下的收敛阶');
legend('u (H^1 误差)', 'ξ (L^2 边界误差)', 'Location', 'best');
grid on;