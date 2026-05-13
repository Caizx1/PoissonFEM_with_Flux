%% plot_meshes_with_boundary_nodes.m
% 绘制 h=0.1 时的内部网格，并独立显示边界节点（黑色点，不连线）
% 线性元：边界网格步长 = 1.25 * h
% 二次元：边界网格步长 = h

clear; clc; close all;

% 几何区域（单位正方形）
geom = Rectg(0, 0, 1, 1);

% 内部网格尺寸
h_in = 0.05;

% 生成内部网格（线性剖分，与基函数阶数无关）
[p0, ~, t0] = initmesh(geom, 'Hmax', h_in);
p = p0';                     % np×2
t = t0(1:3, :)';             % nt×3

% ==================== 线性元情形 ====================
h_bd_linear = 1.25 * h_in;
e_b_nodes_lin = boundaryDivide(geom, h_bd_linear);

figure('Name', '线性元：内部网格与边界节点 (h_bd = 1.25h)', 'Position', [100, 100, 800, 600]);
triplot(t, p(:,1), p(:,2), 'k-', 'LineWidth', 0.5);
hold on;
% 只标记边界节点（黑色圆点），不连线
plot(e_b_nodes_lin(1,:), e_b_nodes_lin(2,:), 'ko', 'MarkerSize', 3, 'LineWidth', 1);
axis equal; grid on;
xlabel('x'); ylabel('y');
% title('线性元 (内部网格 h=0.1, 边界网格步长 1.25h)');
% legend('内部网格', '边界节点 (独立)', 'Location', 'best');

% ==================== 二次元情形 ====================
h_bd_quad = h_in;
e_b_nodes_quad = boundaryDivide(geom, h_bd_quad);

figure('Name', '二次元：内部网格与边界节点 (h_bd = h)', 'Position', [100, 100, 800, 600]);
triplot(t, p(:,1), p(:,2), 'k-', 'LineWidth', 0.5);
hold on;
plot(e_b_nodes_quad(1,:), e_b_nodes_quad(2,:), 'ko', 'MarkerSize', 4, 'LineWidth', 1);
axis equal; grid on;
xlabel('x'); ylabel('y');
% title('二次元 (内部网格 h=0.1, 边界网格步长 h)');
% legend('内部网格', '边界节点 (独立)', 'Location', 'best');