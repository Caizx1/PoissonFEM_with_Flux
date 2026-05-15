function [sol, mesh] = fem_normal_derivative(geom, f, g, h, order, refine_opts)
% 标准有限元求解 Poisson 方程 + 直接法向导数
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
%   sol - 结构体 .u, .xi
%   mesh - 结构体 .p, .t, .e, .xi_location

    if nargin < 5, order = 'linear'; end
    if nargin < 6, refine_opts = struct('use', false); end
    if ~isfield(refine_opts, 'use'), refine_opts.use = false; end
    if ~isfield(refine_opts, 'C'), refine_opts.C = 1.5; end
    if ~isfield(refine_opts, 'max_iter'), refine_opts.max_iter = 10; end

    % 1. 生成内部网格
    [p, t, e] = generate_mesh(geom, h, order, refine_opts);
    % p: np×2, t: nt×3 或 nt×6, e: ne×2 (顶点对)

    % 2. 组装刚度矩阵和载荷向量
    p_for_assembler = p';   % 2×np
    t_for_assembler = t';   % 3×nt 或 6×nt
    A = StiffnessAssembler2D(p_for_assembler, t_for_assembler, @(x,y)1, order);
    F = LoadAssembler2D(p_for_assembler, t_for_assembler, f, order);

    % 3. 施加 Dirichlet 边界条件并求解
    u = apply_dirichlet(A, F, p, t, g, order);

    % 4. 后处理计算边界法向导数及坐标
    if strcmp(order, 'linear')
        [xi, e_out] = project_normal_derivative_p1(u, p, t(:,1:3));
        ne = size(e_out,1);
        xi_loc = zeros(ne,2);
        for i = 1:ne
            v1 = e_out(i,1); v2 = e_out(i,2);
            xi_loc(i,:) = (p(v1,:) + p(v2,:)) / 2;
        end
    else
        % 二次元：节点模式输出
        [xi, e_out, bdNodeId] = project_normal_derivative_p2(u, p, t, 'mode', 'node');
        % xi: ne×3, e_out: ne×2, bdNodeId: ne×3
        ne = size(e_out,1);
        xi_loc = zeros(ne,3,2);
        for i = 1:ne
            nodes = bdNodeId(i,:);
            xi_loc(i,:,:) = p(nodes,:);
        end
    end

    % 5. 构建输出结构体
    sol = struct('u', u, 'xi', xi);
    mesh = struct('p', p, 't', t, 'e', e_out, 'xi_location', xi_loc);
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

function u = apply_dirichlet(A, F, p, t, g, order)
    if strcmp(order, 'linear')
        T = auxstructure(t);
        bdEdges = T.bdEdge;
        bdNodes = unique(bdEdges(:));
    else
        % 二次元：提取所有边界节点（顶点 + 边中点）
        % 1) 获取边界顶点对及所属单元
        T = auxstructure(t(:,1:3));
        bdEdges = T.bdEdge;      % ne×2
        bdEdge2elem = T.bdEdge2elem; % ne×1
        
        bdNodes_vertex = unique(bdEdges(:));
        bdNodes_mid = [];
        % 2) 遍历每条边界边，找到对应的中点索引
        for i = 1:size(bdEdges,1)
            v1 = bdEdges(i,1); v2 = bdEdges(i,2);
            K = bdEdge2elem(i);
            nodesK = t(K,:);
            % 调用子函数获取边上三个节点的顺序（起点、中点、终点）
            [~, edge_loc] = get_edge_nodes_order(v1, v2, nodesK);
            switch edge_loc
                case 1  % 边 2-3，中点索引为 nodesK(5)
                    mid = nodesK(5);
                case 2  % 边 3-1，中点索引为 nodesK(6)
                    mid = nodesK(6);
                case 3  % 边 1-2，中点索引为 nodesK(4)
                    mid = nodesK(4);
                otherwise
                    error('无法确定边号');
            end
            bdNodes_mid = [bdNodes_mid; mid];
        end
        bdNodes = unique([bdNodes_vertex; bdNodes_mid]);
    end
    
    % 修正内部节点的载荷向量
    g_bd = g(p(bdNodes,1), p(bdNodes,2));
    if isscalar(g_bd) && length(bdNodes) > 1
        g_bd = repmat(g_bd, length(bdNodes), 1);
    end
    allNodes = 1:size(p,1);
    inNodes = setdiff(allNodes, bdNodes);
    F(inNodes) = F(inNodes) - A(inNodes, bdNodes) * g_bd;

    nbd = length(bdNodes);
    A(bdNodes,:) = 0;
    A(:,bdNodes) = 0;
    A(bdNodes,bdNodes) = speye(nbd);
    F(bdNodes) = g_bd;
    u = A \ F;
    u = full(u);
end

function [node_order, edge_loc] = get_edge_nodes_order(v1, v2, nodesK)
% 获取二次三角形单元中指定边（由两个顶点 v1,v2 定义）上的三个节点顺序
% 输入：
%   v1, v2   : 边的两个顶点索引
%   nodesK   : 单元节点索引向量 (1×6)
% 输出：
%   node_order : [起点, 中点, 终点] 的全局索引
%   edge_loc   : 边在单元中的局部编号 (1,2,3 对应边 2-3, 3-1, 1-2)
    vK = nodesK(1:3);
    locV1 = find(vK == v1, 1);
    locV2 = find(vK == v2, 1);
    if isempty(locV1) || isempty(locV2)
        error('边界边顶点不属于指定单元');
    end
    if (locV1==1 && locV2==2) || (locV1==2 && locV2==1)
        edge_loc = 3;
        mid_node = nodesK(4);
    elseif (locV1==2 && locV2==3) || (locV1==3 && locV2==2)
        edge_loc = 1;
        mid_node = nodesK(5);
    elseif (locV1==3 && locV2==1) || (locV1==1 && locV2==3)
        edge_loc = 2;
        mid_node = nodesK(6);
    else
        error('无法确定边号');
    end
    % 确保 node_order 顺序为 v1 -> 中点 -> v2
    node_order = [v1, mid_node, v2];
end