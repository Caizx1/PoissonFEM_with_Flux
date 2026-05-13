%% test_ratio_convergence2.m
% 使用 H1_error 和 L2_boundary_error 函数，研究不同边界网格步长比下二次内部元+线性边界元的收敛阶
% 输出每个 ratio 下 u_H1 误差和 xi_L2 误差的收敛阶，并绘制收敛阶随 ratio 变化的曲线

clear; clc; close all;

% ===================== 1. 问题定义 =====================

% geom = Rectg(0, 0, 1, 1);                % 几何
% u_exact = @(x,y) sin(pi*x) .* sin(pi*y); % 精确解
% grad_u_exact = @(x,y) [pi*cos(pi*x).*sin(pi*y), pi*sin(pi*x).*cos(pi*y)]; % 精确梯度
% f = @(x,y) 2*pi^2 * sin(pi*x) .* sin(pi*y); % 右端项
% g = @(x,y) 0;                            % 齐次 Dirichlet
% 
% xi_exact_fun = @(x,y) ...                % 精确边界法向导数 ξ = -∂u/∂n
%     (abs(x)<1e-12) * (-pi * sin(pi*y)) + ...
%     (abs(x-1)<1e-12) * (-pi * sin(pi*y)) + ...
%     (abs(y)<1e-12) * (-pi * sin(pi*x)) + ...
%     (abs(y-1)<1e-12) * (-pi * sin(pi*x));

% ===================== 非齐次 Dirichlet 算例（纯指数-三角） =====================
% 精确解：u(x,y) = exp(x) * sin(pi*y) + exp(y) * sin(pi*x)
% 该解在单位正方形边界上非零，光滑且不含多项式，右端项 f = (pi^2-1) * u

% 几何：单位正方形
geom = Rectg(0, 0, 1, 1);

% 精确解
u_exact = @(x,y) exp(x) .* sin(pi*y) + exp(y) .* sin(pi*x);

% 精确梯度
grad_u_exact = @(x,y) [exp(x).*sin(pi*y) + pi*exp(y).*cos(pi*x), ...
                       pi*exp(x).*cos(pi*y) + exp(y).*sin(pi*x)];

% 右端项 f = -Δu = (pi^2 - 1) * u
f = @(x,y) (pi^2 - 1) * (exp(x).*sin(pi*y) + exp(y).*sin(pi*x));

% Dirichlet 边界条件
g = @(x,y) u_exact(x,y);

xi_exact_fun = @(x,y) ...
    (abs(y)   < 1e-12) * ( -pi*exp(x) - sin(pi*x) ) + ...          % 下边界
    (abs(y-1) < 1e-12) * ( -pi*exp(x) + exp(1)*sin(pi*x) ) + ...   % 上边界
    (abs(x)   < 1e-12) * ( -sin(pi*y) - pi*exp(y) ) + ...           % 左边界
    (abs(x-1) < 1e-12) * ( exp(1)*sin(pi*y) - pi*exp(y) );          % 右边界

% % ===================== 非齐次 Dirichlet 算例（法向导数连续） =====================
% % 精确解：u(x,y) = sin(π x) + sin(π y)
% % 该解在单位正方形边界上非零，且法向导数在角点处连续
% 
% % 几何：单位正方形
% geom = Rectg(0, 0, 1, 1);
% % 精确解
% u_exact = @(x,y) sin(pi*x) + sin(pi*y);
% % 精确梯度
% grad_u_exact = @(x,y) [pi*cos(pi*x), pi*cos(pi*y)];
% % 右端项 f = -Δu = π^2 sin(π x) + π^2 sin(π y)
% f = @(x,y) pi^2 * (sin(pi*x) + sin(pi*y));
% % Dirichlet 边界条件
% g = @(x,y) u_exact(x,y);
% xi_exact_fun = @(x,y) -pi * ones(size(x));


% ===================== 2. 参数设置 =====================
h_list = [0.2, 0.1, 0.05, 0.025];       % 内部网格尺寸（粗到细）
ratio_list = [0.5, 0.75, 1, 1.25, 1.5, 2]; % 边界/内部步长比

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
        [u, xi, p, t, e, e_b_nodes, e_boundary] = ...
            primal_mixed_solver2D(geom, f, g, h_in, h_bd, 'quadratic', refine_opts);
        
        % 调整数据格式以适应误差函数
        p = p';               % 2×np → np×2
        t = t';               % 6×nt → nt×6
        
        % 计算 u 的 H1 误差（需要梯度函数）
        err_u_val = H1_error(p, t, u, u_exact, grad_u_exact, 'quadratic');
        err_u = [err_u, err_u_val];
        
        % 计算 xi 的 L2 边界误差（使用独立边界网格信息）
        % 注意：e_boundary 是 2×ne，e_b_nodes 是 2×nb
        err_xi_val = L2_boundary_error(e_b_nodes, e_boundary, xi, xi_exact_fun);
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
title('二次内部元 + 线性边界元：不同边界网格步长比下的收敛阶');
legend('u (H^1 误差)', 'ξ (L^2 边界误差)', 'Location', 'best');
grid on;

% 可选：绘制误差随 h 变化的详细曲线（最后一个比值）
% 这里略去，可按需添加