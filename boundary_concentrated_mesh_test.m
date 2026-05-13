%% 测试 boundary_concentrated_mesh 函数
% 需要确保 bisect 和 auxstructure 在 MATLAB 路径中（ifem 包）

clear; clc; close all;

%% 1. 定义几何（矩形区域 [0,1]×[0,1]）
r = Rectg(0, 0, 1, 1);   % 几何矩阵

%% 2. 生成初始粗网格
hInit = 0.2;                % 初始最大单元直径
[p0, e0, elem0] = initmesh(r, 'Hmax', hInit);
p0 = p0';                   % 转换为 np×2
elem0 = elem0(1:3,:)';      % 转换为 nt×3

% 计算初始内部最大直径（作为 h 输入）
h = hInit;

%% 3. 边界集中加密
C = 1.5;                    % 判据常数
max_iter = 15;              % 最大迭代次数
[p, e, elem] = boundary_concentrated_mesh(p0, elem0, h, r, C, max_iter);

%% 4. 计算加密后网格的单元直径和到边界的距离
hT = elem_diam(p, elem);
rho = dist2bd(p, elem, boundary_segment(r));

%% 5. 可视化结果
% 5.1 显示加密后的网格
figure('Name', '加密后网格', 'Position', [100, 100, 800, 600]);
triplot(elem, p(:,1), p(:,2), 'k-', 'LineWidth', 0.5);
axis equal; hold on;
% 高亮边界边
plot([p(e(1,:),1); p(e(2,:),1)], [p(e(1,:),2); p(e(2,:),2)], 'r-', 'LineWidth', 2);
title('边界集中网格（红色为边界）');
xlabel('x'); ylabel('y');

% 5.2 单元直径与到边界距离的散点图（对数坐标）
figure('Name', '网格尺寸分布', 'Position', [100, 100, 800, 600]);
loglog(rho, hT, 'b.', 'MarkerSize', 8);
hold on;
% 理论曲线
rho_plot = logspace(log10(min(rho(rho>0))), log10(max(rho)), 100);
plot(rho_plot, C * h * sqrt(rho_plot), 'r--', 'LineWidth', 2);
plot(rho_plot, C * h^2 * ones(size(rho_plot)), 'k--', 'LineWidth', 2);
xlabel('单元到边界距离 ρ');
ylabel('单元直径 h_T');
title('网格尺寸分布与理论判据');
legend('实际单元', '内部判据 h_T ~ C h √ρ', '边界判据 h_T ~ C h²', 'Location', 'best');
grid on;

% 5.3 单元到边界距离的云图（可选）
figure('Name', '单元到边界距离云图', 'Position', [100, 100, 800, 600]);
patch('Faces', elem, 'Vertices', p, 'FaceVertexCData', rho, ...
      'FaceColor', 'flat', 'EdgeColor', 'none');
view(2); axis equal; colorbar;
xlabel('x'); ylabel('y'); title('单元到边界距离 ρ');

%% 6. 统计信息
fprintf('初始网格: 节点数 = %d, 单元数 = %d\n', size(p0,1), size(elem0,1));
fprintf('加密后网格: 节点数 = %d, 单元数 = %d\n', size(p,1), size(elem,1));
fprintf('内部最大直径 h = %g\n', h);
fprintf('边界最小单元直径 ≈ %g\n', min(hT(rho < 1e-8)));
fprintf('边界最大单元直径 ≈ %g\n', max(hT(rho < 1e-8)));
fprintf('理论边界目标尺寸 C*h^2 = %g\n', C * h^2);