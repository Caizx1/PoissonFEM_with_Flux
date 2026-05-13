% test_boundary_normal_convergence.m
% 二次有限元边界法向导数收敛性测试

clear; clc; close all;

%% 问题设置
geom = Rectg(0,0,1,1);
u_exact = @(x,y) sin(pi*x) .* sin(pi*y);
a = @(x,y) 1;
f = @(x,y) 2*pi^2 * sin(pi*x) .* sin(pi*y);
g = u_exact;

% 网格尺寸列表
h_list = [0.2, 0.1, 0.05, 0.025,0.0125];
n_h = length(h_list);
L2_errs = zeros(n_h, 1);
h_actual = zeros(n_h, 1);  % 实际平均边长

%% 循环不同网格尺寸
fprintf('\n 二次基函数计算投影法向导数');
for i_h = 1:n_h
    h_target = h_list(i_h);
    fprintf('\n========== h = %.4f ==========\n', h_target);

    % --- 生成初始线性网格 ---
    [pl, ~, tl] = initmesh(geom, 'Hmax', h_target);
    p_lin = pl';                  % np×2
    t_lin = tl(1:3,:)';           % nt×3

    % --- 转换为二次网格（6节点三角形）---
    [p_quad, t_quad] = linear2quadMesh(p_lin, t_lin);
    p = p_quad';                  % 2×np_quad
    t = t_quad';                  % 6×nt

    % 记录实际网格尺寸（平均边界边长）
    Ttmp = auxstructure(t(1:3,:)');
    bdEdge_tmp = Ttmp.bdEdge;
    edge_lengths = zeros(size(bdEdge_tmp,1),1);
    for ie = 1:size(bdEdge_tmp,1)
        v1 = bdEdge_tmp(ie,1); v2 = bdEdge_tmp(ie,2);
        edge_lengths(ie) = norm(p(:,v1)-p(:,v2));
    end
    h_actual(i_h) = mean(edge_lengths);

    % --- 有限元组装 ---
    A = StiffnessAssembler2D(p, t, a, 'quadratic');
    b = LoadAssembler2D(p, t, f, 'quadratic');

    % --- 提取所有 Dirichlet 边界节点（含中点）---
    Tstruct = auxstructure(t(1:3,:)');
    bdEdge = Tstruct.bdEdge;
    bdEdge2elem = Tstruct.bdEdge2elem;

    bdNodes_vertex = unique(bdEdge(:));
    bdNodes_mid = [];
    for i = 1:size(bdEdge,1)
        v1 = bdEdge(i,1); v2 = bdEdge(i,2);
        K = bdEdge2elem(i);
        nodesK = t(:,K)';
        [~, edge_loc] = get_edge_nodes_order(v1, v2, nodesK);
        switch edge_loc
            case 1  % 边 2-3，中点5
                mid = nodesK(5);
            case 2  % 边 3-1，中点6
                mid = nodesK(6);
            case 3  % 边 1-2，中点4
                mid = nodesK(4);
            otherwise
                error('无法确定边号');
        end
        bdNodes_mid = [bdNodes_mid; mid];
    end
    bdNodes = unique([bdNodes_vertex; bdNodes_mid]);
    nbd = length(bdNodes);

    % 施加 Dirichlet 边界条件
    A(bdNodes,:) = 0;
    A(:,bdNodes) = 0;
    A(bdNodes,bdNodes) = speye(nbd);
    b(bdNodes) = g(p(1,bdNodes), p(2,bdNodes));

    u = A \ b;

    % --- 计算边界法向导数（高斯点模式）---
    gauss = project_normal_derivative_p2(u, p', t', 'mode', 'gauss', 'nGauss', 3);

    % --- 计算 L2 误差（利用高斯点所在边的法向）---
    nPts = size(gauss.x, 1);
    exact_vals = zeros(nPts, 1);
    for i = 1:nPts
        eid = gauss.edgeIdx(i);
        v1 = bdEdge(eid,1); v2 = bdEdge(eid,2);
        K = bdEdge2elem(eid);
        xyK = p(:,t(1:3,K))';
        n_outer = outer_normal(xyK, p(:,[v1,v2])');
        x = gauss.x(i,1); y = gauss.x(i,2);
        grad_ex = [pi*cos(pi*x)*sin(pi*y); pi*sin(pi*x)*cos(pi*y)];
        exact_vals(i) = dot(grad_ex, n_outer);
    end

    L2_err_sq = sum(gauss.weight .* (gauss.dudn - exact_vals).^2);
    L2_errs(i_h) = sqrt(L2_err_sq);
    fprintf('边界法向导数 L2 误差 = %e\n', L2_errs(i_h));
end

%% 计算收敛阶
fprintf('\n========== 收敛性分析 ==========\n');
fprintf('  h_avg        L2_error        order\n');
for i_h = 1:n_h
    if i_h == 1
        fprintf('%.4f      %.6e      ---\n', h_actual(i_h), L2_errs(i_h));
    else
        order = log2(L2_errs(i_h-1) / L2_errs(i_h));
        fprintf('%.4f      %.6e      %.2f\n', h_actual(i_h), L2_errs(i_h), order);
    end
end
order_avg = polyfit(log(h_list),log(L2_errs),1);
order_avg = order_avg(1);
fprintf('\n 收敛阶 = %.2f\n',order_avg);

% 绘制误差图
figure;
loglog(h_actual, L2_errs, 'bo-', 'LineWidth', 1.5);
hold on;
loglog(h_actual, h_actual.^2 * (L2_errs(1)/h_actual(1)^2), 'r--', 'DisplayName', 'O(h^2)');
xlabel('平均边界边长 h');
ylabel('边界法向导数 L2 误差');
legend('数值误差', '二阶参考');
grid on;
title('二次元边界法向导数收敛性');

%% 子函数定义（与 project_normal_derivative_p2.m 中保持一致）
function [node_order, edge_loc] = get_edge_nodes_order(v1, v2, nodesK)
    vK = nodesK(1:3);
    locV1 = find(vK == v1, 1);
    locV2 = find(vK == v2, 1);
    if isempty(locV1) || isempty(locV2)
        error('边界边顶点不属于指定单元');
    end
    if (locV1==1 && locV2==2) || (locV1==2 && locV2==1)
        edge_loc = 3; mid_node = nodesK(4);
    elseif (locV1==2 && locV2==3) || (locV1==3 && locV2==2)
        edge_loc = 1; mid_node = nodesK(5);
    elseif (locV1==3 && locV2==1) || (locV1==1 && locV2==3)
        edge_loc = 2; mid_node = nodesK(6);
    else
        error('无法确定边号');
    end
    node_order = [v1, mid_node, v2];
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