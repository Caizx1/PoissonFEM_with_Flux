function [sol, mesh] = fem_variational_derivative(geom, f, g, h, order, refine_opts)
% 变分法计算法向导数：通过残差和边界质量矩阵求解 λ = ∂u/∂n
% 求解 Poisson 方程 -Δu = f，Dirichlet 边界 u = g
%
% 输入：
%   geom       - 几何描述
%   f          - 右端项函数句柄
%   g          - Dirichlet 边界条件函数句柄
%   h          - 网格尺寸
%   order      - 'linear' 或 'quadratic'，默认 'linear'
%   refine_opts- 边界加密选项（可选，结构体）
%
% 输出：
%   sol - 结构体 .u (节点解), .lambda (边界节点处的法向导数)
%   mesh - 结构体 .p, .t, .e (边界顶点对), .lambda_location (边界节点坐标)

    if nargin < 5, order = 'linear'; end
    if nargin < 6, refine_opts = struct('use', false); end
    if ~isfield(refine_opts, 'use'), refine_opts.use = false; end
    if ~isfield(refine_opts, 'C'), refine_opts.C = 1.5; end
    if ~isfield(refine_opts, 'max_iter'), refine_opts.max_iter = 10; end

    % 1. 生成内部网格
    [p, t, e] = generate_mesh(geom, h, order, refine_opts);
    % p: np×2, t: nt×3 或 nt×6, e: ne×2 (边界顶点对)

    % 2. 组装原始刚度矩阵和载荷向量（不施加边界条件）
    p_for_assembler = p';   % 2×np
    t_for_assembler = t';   % 3×nt 或 6×nt
    A_orig = StiffnessAssembler2D(p_for_assembler, t_for_assembler, @(x,y)1, order);
    F_orig = LoadAssembler2D(p_for_assembler, t_for_assembler, f, order);

    % 3. 施加 Dirichlet 边界条件（正确修改右端项）并求解 u
    [A_mod, F_mod, bdNodes] = apply_dirichlet(A_orig, F_orig, p, t, g, order);
    u = A_mod \ F_mod;
    u = full(u);

    % 4. 计算残差 r = A_orig * u - F_orig
    r = A_orig * u - F_orig;

    % 5. 提取边界节点上的残差 r_b，组装边界质量矩阵 M_b，求解 lambda
    if strcmp(order, 'linear')
        % 线性元：边界节点即边界边的顶点（注意每个顶点可能被多条边共享，但质量矩阵组装按边进行，最后求解的是节点值）
        bdNodes = unique(e(:));  % 边界节点列表
        n_bd = length(bdNodes);
        r_b = r(bdNodes);
        M_b = boundary_mass_matrix_linear(p, e, bdNodes);
        lambda = M_b \ r_b;
        % 输出位置：边界节点坐标
        lambda_loc = p(bdNodes, :);
    else
        % 二次元：边界节点包括顶点和边中点，需要正确收集所有边界节点
        % 首先获取边界边信息（顶点对）及其中点
        T = auxstructure(t(:,1:3));
        bdEdges = T.bdEdge;          % ne×2
        bdEdge2elem = T.bdEdge2elem; % ne×1
        ne = size(bdEdges,1);
        % 收集所有边界节点（顶点+中点），并记录每条边对应的三个节点
        bdNodes_vertex = unique(bdEdges(:));
        bdNodes_mid = zeros(ne,1);
        edge_nodes = zeros(ne,3);    % 每行 [v1, mid, v2]
        for i = 1:ne
            v1 = bdEdges(i,1); v2 = bdEdges(i,2);
            K = bdEdge2elem(i);
            nodesK = t(K,:);
            [node_order, ~] = get_edge_nodes_order(v1, v2, nodesK);
            edge_nodes(i,:) = node_order;  % [v1, mid, v2]
            bdNodes_mid(i) = node_order(2);
        end
        bdNodes = unique([bdNodes_vertex; bdNodes_mid]);
        % 建立全局节点到局部编号的映射
        [~, bdNodeId] = ismember(bdNodes, bdNodes);
        % 但组装质量矩阵需要按边进行，使用 edge_nodes
        n_bd = length(bdNodes);
        r_b = r(bdNodes);
        M_b = boundary_mass_matrix_quadratic(p, edge_nodes, bdNodes);
        lambda = M_b \ r_b;
        % 输出位置：边界节点的坐标（顺序同 bdNodes）
        lambda_loc = p(bdNodes, :);
    end

    % 6. 构建输出结构体
    sol = struct('u', u, 'lambda', lambda);
    mesh = struct('p', p, 't', t, 'e', e, 'lambda_location', lambda_loc);
end

% ==================== 内部子函数 ====================

function [p, t, e] = generate_mesh(geom, h, order, refine_opts)
    [p0, ~, t0] = initmesh(geom, 'Hmax', h);
    p0 = p0';   % np0×2
    t0 = t0(1:3,:)';  % nt0×3
    if refine_opts.use
        [p, e, t] = boundary_concentrated_mesh(p0, t0, h, geom, ...
            refine_opts.C, refine_opts.max_iter);
        e = e';
    else
        p = p0;
        t = t0;
        T = auxstructure(t);
        e = T.bdEdge;
    end
    if strcmp(order, 'quadratic')
        [p, t] = linear2quadMesh(p, t);
        % 重新提取边界边（顶点对）
        T = auxstructure(t(:,1:3));
        e = T.bdEdge;
    end
end

function [A_mod, F_mod, bdNodes] = apply_dirichlet(A, F, p, t, g, order)
% 施加 Dirichlet 边界条件，同时修改右端项以保持内部方程的正确性
    if strcmp(order, 'linear')
        T = auxstructure(t);
        bdEdges = T.bdEdge;
        bdNodes = unique(bdEdges(:));
    else
        % 二次元：提取所有边界节点（顶点 + 边中点）
        T = auxstructure(t(:,1:3));
        bdEdges = T.bdEdge;
        bdEdge2elem = T.bdEdge2elem;
        bdNodes_vertex = unique(bdEdges(:));
        bdNodes_mid = [];
        for i = 1:size(bdEdges,1)
            v1 = bdEdges(i,1); v2 = bdEdges(i,2);
            K = bdEdge2elem(i);
            nodesK = t(K,:);
            [~, edge_loc] = get_edge_nodes_order(v1, v2, nodesK);
            switch edge_loc
                case 1, mid = nodesK(5);
                case 2, mid = nodesK(6);
                case 3, mid = nodesK(4);
            end
            bdNodes_mid = [bdNodes_mid; mid];
        end
        bdNodes = unique([bdNodes_vertex; bdNodes_mid]);
    end

    g_bd = g(p(bdNodes,1), p(bdNodes,2));
    if isscalar(g_bd) && length(bdNodes) > 1
        g_bd = repmat(g_bd, length(bdNodes), 1);
    end
    % 内部节点修正
    allNodes = 1:size(p,1);
    inNodes = setdiff(allNodes, bdNodes);
    F(inNodes) = F(inNodes) - A(inNodes, bdNodes) * g_bd;

    A_mod = A;
    F_mod = F;
    A_mod(bdNodes,:) = 0;
    A_mod(:,bdNodes) = 0;
    A_mod(bdNodes,bdNodes) = speye(length(bdNodes));
    F_mod(bdNodes) = g_bd;
end

function [node_order, edge_loc] = get_edge_nodes_order(v1, v2, nodesK)
    vK = nodesK(1:3);
    locV1 = find(vK == v1, 1);
    locV2 = find(vK == v2, 1);
    if (locV1==1 && locV2==2) || (locV1==2 && locV2==1)
        edge_loc = 3; mid_node = nodesK(4);
    elseif (locV1==2 && locV2==3) || (locV1==3 && locV2==2)
        edge_loc = 1; mid_node = nodesK(5);
    elseif (locV1==3 && locV2==1) || (locV1==1 && locV2==3)
        edge_loc = 2; mid_node = nodesK(6);
    else
        error('边顶点不属于单元');
    end
    node_order = [v1, mid_node, v2];
end

function M = boundary_mass_matrix_linear(p, e, bdNodes)
% 线性元边界质量矩阵（仅考虑边界边）
% p: np×2, e: ne×2 (边界顶点对), bdNodes: 边界节点列表
    ne = size(e,1);
    M = sparse(length(bdNodes), length(bdNodes));
    node2idx = zeros(size(p,1),1);
    node2idx(bdNodes) = 1:length(bdNodes);
    for i = 1:ne
        n1 = e(i,1); n2 = e(i,2);
        idx1 = node2idx(n1); idx2 = node2idx(n2);
        if idx1==0 || idx2==0, error('边界节点不在bdNodes中'); end
        p1 = p(n1,:); p2 = p(n2,:);
        L = norm(p2 - p1);
        Me = [2,1;1,2] * L / 6;
        M([idx1,idx2],[idx1,idx2]) = M([idx1,idx2],[idx1,idx2]) + Me;
    end
end

function M = boundary_mass_matrix_quadratic(p, edge_nodes, bdNodes)
% 二次元边界质量矩阵
% edge_nodes: ne×3，每行 [v1, mid, v2]
% bdNodes: 所有边界节点的全局索引列表（顶点+中点）
    ne = size(edge_nodes,1);
    M = sparse(length(bdNodes), length(bdNodes));
    node2idx = zeros(size(p,1),1);
    node2idx(bdNodes) = 1:length(bdNodes);
    for i = 1:ne
        v1 = edge_nodes(i,1); mid = edge_nodes(i,2); v2 = edge_nodes(i,3);
        idx1 = node2idx(v1); idx2 = node2idx(v2); idx3 = node2idx(mid);
        if any([idx1,idx2,idx3]==0), error('节点不在bdNodes中'); end
        p1 = p(v1,:); p2 = p(v2,:);
        L = norm(p2 - p1);
        Me = [4, -1, 2; -1, 4, 2; 2, 2, 16] * L / 30;
        M([idx1,idx2,idx3],[idx1,idx2,idx3]) = M([idx1,idx2,idx3],[idx1,idx2,idx3]) + Me;
    end
end