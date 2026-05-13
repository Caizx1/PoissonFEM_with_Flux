% test_linear_convergence_clean.m
% 线性元边界法向导数收敛性测试（直接投影法计算）

clear; clc; close all;

%% 问题设置
geom = Rectg(0,0,1,1);
u_exact = @(x,y) sin(pi*x) .* sin(pi*y);
a = @(x,y) 1;
f = @(x,y) 2*pi^2 * sin(pi*x) .* sin(pi*y);
g = u_exact;

h_list = [0.2, 0.1, 0.05, 0.025, 0.0125];
n_h = length(h_list);
L2_errs = zeros(n_h, 1);
h_actual = zeros(n_h, 1);

fprintf('\n 线性基函数计算投影法向导数')
for i_h = 1:n_h
    h_target = h_list(i_h);
    fprintf('\n========== h = %.4f ==========\n', h_target);

    % --- 生成网格 ---
    [pl, ~, tl] = initmesh(geom, 'Hmax', h_target);
    p = pl';                     % np×2
    t = tl(1:3,:)';              % nt×3

    % 统一为逆时针定向
    for k = 1:size(t,1)
        xy = p(t(k,:),:);
        if polyarea(xy(:,1), xy(:,2)) < 0
            t(k,[2,3]) = t(k,[3,2]);
        end
    end

    % 记录平均边界边长
    Ttmp = auxstructure(t);
    bdEdge_tmp = Ttmp.bdEdge;
    edge_lengths = zeros(size(bdEdge_tmp,1),1);
    for ie = 1:size(bdEdge_tmp,1)
        v1 = bdEdge_tmp(ie,1); v2 = bdEdge_tmp(ie,2);
        edge_lengths(ie) = norm(p(v1,:)-p(v2,:));
    end
    h_actual(i_h) = mean(edge_lengths);

    % --- 有限元组装 ---
    A = StiffnessAssembler2D(p', t', a, 'linear');
    b = LoadAssembler2D(p', t', f, 'linear');

    % 提取边界节点并施加 Dirichlet 条件
    Tstruct = auxstructure(t);
    bdEdge = Tstruct.bdEdge;
    bdEdge2elem = Tstruct.bdEdge2elem;
    bdNodes = unique(bdEdge(:));
    nbd = length(bdNodes);

    A(bdNodes,:) = 0;
    A(:,bdNodes) = 0;
    A(bdNodes,bdNodes) = speye(nbd);
    b(bdNodes) = g(p(bdNodes,1), p(bdNodes,2));

    u = A \ b;

    % --- 计算边界法向导数 ---
    ne = size(bdEdge,1);
    xi = zeros(ne,1);
    for ie = 1:ne
        v1 = bdEdge(ie,1); v2 = bdEdge(ie,2);
        K = bdEdge2elem(ie);
        nodes = t(K,:);
        xy = p(nodes,:);
        [~, bvec, cvec] = HatGradients(xy(:,1), xy(:,2));
        grad_u = [bvec, cvec]' * u(nodes);
        n_outer = outer_normal(xy, p([v1;v2],:));
        xi(ie) = dot(grad_u, n_outer);
    end

    % --- L2 误差计算 ---
    err_sq = zeros(ne,1);
    for ie = 1:ne
        v1 = bdEdge(ie,1); v2 = bdEdge(ie,2);
        xm = 0.5 * (p(v1,1) + p(v2,1));
        ym = 0.5 * (p(v1,2) + p(v2,2));
        if xm < 1e-8
            exact_val = -pi * sin(pi*ym);
        elseif xm > 1-1e-8
            exact_val = -pi * sin(pi*ym);
        elseif ym < 1e-8
            exact_val = -pi * sin(pi*xm);
        elseif ym > 1-1e-8
            exact_val = -pi * sin(pi*xm);
        else
            error('边中点不在边界上');
        end
        he = norm(p(v2,:) - p(v1,:));
        err_sq(ie) = he * (xi(ie) - exact_val)^2;
    end
    L2_err = sqrt(sum(err_sq));
    L2_errs(i_h) = L2_err;
    fprintf('边界法向导数 L2 误差 = %e\n', L2_err);
end

%% 收敛性分析
fprintf('\n========== 收敛性分析 ==========\n');
fprintf('  h_avg        L2_error        order\n');
for i_h = 1:n_h
    if i_h == 1
        fprintf('%.4f      %.6e      ---\n', h_actual(i_h), L2_errs(i_h));
    else
        order = log2(L2_errs(i_h-1) / L2_errs(i_h)) / log2(h_actual(i_h-1) / h_actual(i_h));
        fprintf('%.4f      %.6e      %.2f\n', h_actual(i_h), L2_errs(i_h), order);
    end
end
order_avg = polyfit(log(h_list),log(L2_errs),1);
order_avg = order_avg(1);
fprintf('\n 收敛阶 = %.2f\n',order_avg);

% 绘图
figure;
loglog(h_actual, L2_errs, 'bo-', 'LineWidth', 1.5);
hold on;
loglog(h_actual, h_actual * (L2_errs(1)/h_actual(1)), 'r--', 'DisplayName', 'O(h)');
xlabel('平均边界边长 h');
ylabel('边界法向导数 L2 误差');
legend('数值误差', '一阶参考', 'Location', 'best');
grid on;
title('线性元边界法向导数收敛性');

%% 辅助函数 ----------------------------------------------------------------
function [area, b, c] = HatGradients(x, y)
    x = x(:); y = y(:);
    area = 0.5 * abs( x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2)) );
    b = zeros(3,1);
    c = zeros(3,1);
    b(1) = (y(2) - y(3)) / (2*area);
    b(2) = (y(3) - y(1)) / (2*area);
    b(3) = (y(1) - y(2)) / (2*area);
    c(1) = (x(3) - x(2)) / (2*area);
    c(2) = (x(1) - x(3)) / (2*area);
    c(3) = (x(2) - x(1)) / (2*area);
end

function n_out = outer_normal(triangle_vertices, edge_endpoints)
    A = triangle_vertices(1,:); B = triangle_vertices(2,:); C = triangle_vertices(3,:);
    P1 = edge_endpoints(1,:); P2 = edge_endpoints(2,:);
    verts = [A; B; C];
    tol = 1e-12;
    idx = true(3,1);
    for i = 1:3
        if norm(verts(i,:) - P1) < tol || norm(verts(i,:) - P2) < tol
            idx(i) = false;
        end
    end
    third = verts(idx, :);
    d = P2 - P1;
    L = norm(d);
    if L < tol, error('边长度为零'); end
    d_unit = d / L;
    n = [-d_unit(2), d_unit(1)];
    v = P1 - third;
    if dot(n, v) > 0
        n_out = n(:);
    else
        n_out = -n(:);
    end
end