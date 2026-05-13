clear;

r = Rectg(0,0,1,1);
h = 0.08;
[p,~,t] = initmesh(r,"Hmax",h);
p = p';t = t(1:3,:)';
[p,e,t] = boundary_concentrated_mesh(p,t,h,r,1.5,10);
p = p';t = t';
e_b = boundaryDivide(r,0.04 );

A = StiffnessAssembler2D(p,t,@(x,y)1);
F = LoadAssembler2D(p,t,@f);
%[B,G] = BoundaryCoupling2D(p,e,@g);
[B,G] = BoundaryCoupling2D_independent(p,e,e_b,@g);

np = size(p,2);ne_b = size(e_b,2);
u = [A B;B' sparse(ne_b,ne_b)] \ [F;G];

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

% 独立边界网格的边中点
edge_midpoints = (e_b(:, 1:end-1) + e_b(:, 2:end)) / 2;
last_mid = (e_b(:, end) + e_b(:, 1)) / 2;
edge_midpoints = [edge_midpoints, last_mid];

% 绘制边界通量
figure;
scatter(edge_midpoints(1,:), edge_midpoints(2,:), 50, xi, 'filled');
colorbar; xlabel('x'); ylabel('y'); title('Numerical flux ξ_h on boundary');
axis equal;

% 按边序号对比
s = 1:ne_b;
figure;
plot(s, xi, 'b-o', 'LineWidth', 1.5); hold on;

% 计算精确ξ在边中点处的值
xi_exact = zeros(size(s));
for i = 1:length(s)
    m = edge_midpoints(:,i);
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

%=========================

function z = f(x,y)
    z = 2 * pi ^ 2 * sin(pi * x) * sin(pi * y);
end

function z = g(x,y)
    z = 0;
end