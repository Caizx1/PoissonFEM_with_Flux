clear;

r = Rectg(0,0,1,1);
[p,e,t] = initmesh(r,"Hmax",0.01);


A = StiffnessAssembler2D(p,t,@(x,y)1);
F = LoadAssembler2D(p,t,@f);
[B,G] = BoundaryCoupling2D0(p,e,@g);

np = size(p,2);ne = size(e,2);
u = [A B;B' sparse(ne,ne)] \ [F;G];

xi = -u(np + 1:end); xi = full(xi);
u = u(1:np);    u = full(u);



figure(1);
trisurf(t', p(1,:), p(2,:), u, 'EdgeColor', 'none');
colorbar;
xlabel('x'); ylabel('y'); title('Numerical solution u_h');
view(2); axis equal; shading interp;

% 精确解函数
u_exact = @(x,y) sin(pi*x) .* sin(pi*y);

% 在网格节点上计算精确值
U_exact = u_exact(p(1,:), p(2,:))';

figure;
trisurf(t', p(1,:), p(2,:), U_exact, 'EdgeColor', 'none');
colorbar;
xlabel('x'); ylabel('y'); title('Exact solution u');
view(2); axis equal; shading interp;

error_u = U_exact - u;

figure;
trisurf(t', p(1,:), p(2,:), error_u, 'EdgeColor', 'none');
colorbar;
xlabel('x'); ylabel('y'); title('Error u - u_h');
view(2); axis equal; shading interp;


% 计算B的零空间（一维）
eta = null(full(B));          % eta是80×1的向量，满足B*eta≈0

% 从xi中减去零空间分量
xi = xi - (eta' * xi) / (eta' * eta) * eta;


% 绘制每条边界边的中点坐标，用颜色表示ξ_h
midpoints = (p(:, e(1,:)) + p(:, e(2,:))) / 2;

figure;
scatter(midpoints(1,:), midpoints(2,:), 50, xi, 'filled');
colorbar;
xlabel('x'); ylabel('y'); title('Numerical flux ξ_h on boundary');
axis equal;

% 若要沿边界参数s绘制ξ_h与精确值对比，需将边界分段排序
% 简单示例：按边序号画阶梯图（假设边按逆时针顺序）
s = 1:length(xi);  % 边序号
figure;
plot(s, xi, 'b-o', 'LineWidth', 1.5); hold on;
% 计算精确ξ在边中点处的值
xi_exact = zeros(size(s));
for i = 1:length(s)
    m = midpoints(:,i);
    if abs(m(1)) < 1e-12   % 左边
        xi_exact(i) = -pi * sin(pi * m(2));
    elseif abs(m(1)-1) < 1e-12 % 右边
        xi_exact(i) = -pi * sin(pi * m(2));
    elseif abs(m(2)) < 1e-12   % 下边
        xi_exact(i) = -pi * sin(pi * m(1));
    else  % 上边
        xi_exact(i) = -pi * sin(pi * m(1));
    end
end
plot(s, xi_exact, 'r--x', 'LineWidth', 1.5);
xlabel('Boundary edge index'); ylabel('ξ');
legend('ξ_h', 'ξ_{exact}'); title('Boundary flux comparison');



%===========================
% 计算B的零空间（一维）
% eta = null(full(B));          % eta是80×1的向量，满足B*eta≈0

% 从xi中减去零空间分量
% xi_corrected = xi - (eta' * xi) / (eta' * eta) * eta;

% 用xi_corrected替换原来的xi进行误差计算
% xi_error = xi_corrected - xi_exact';
% mean_abs_err = mean(abs(xi_error));
% max_abs_err = max(abs(xi_error));
% fprintf('修正后：平均绝对误差 = %e，最大绝对误差 = %e\n', mean_abs_err, max_abs_err);

%=========================

function z = f(x,y)
    z = 2 * pi ^ 2 * sin(pi * x) * sin(pi * y);
end

function y = g(x)
    y = 0;
end