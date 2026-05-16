function varargout = project_normal_derivative_p2(uh, p, t, options)
% 计算二次三角形单元边界上的法向导数，支持节点输出或高斯点输出
% 输入：
%   uh      : 节点解向量 np×1
%   p       : 节点坐标 np×2
%   t       : 三角形单元 nt×6（顶点+边中点，按对边顺序）
%   options : 可选参数结构体，包含字段：
%             .mode    : 'node' (默认) 或 'gauss'
%             .nGauss  : 高斯点数，默认 3（若 mode='gauss'）
% 输出：
%   若 mode='node':
%       xi       : ne×3，每条边三个节点处的法向导数
%       bdEdge   : ne×2，边界边顶点编号
%       bdNodeId : ne×3，边界边上三个节点的全局编号
%   若 mode='gauss':
%       gauss    : 结构体，包含字段：
%           .x      : 高斯点物理坐标 (nTot×2)
%           .dudn   : 高斯点法向导数值 (nTot×1)
%           .weight : 高斯点积分权重 (nTot×1)（已乘边长 Jacobian）
%           .edgeIdx: 高斯点所属边界边索引 (nTot×1)

    arguments
        uh (:,1) double
        p (:,2) double
        t (:,6) double
        options.mode char {mustBeMember(options.mode,{'node','gauss'})} = 'node'
        options.nGauss (1,1) double {mustBeInteger, mustBePositive} = 3
    end

    % 使用顶点子网格提取边界边信息
    T = auxstructure(t(:,1:3));
    bdEdge = T.bdEdge;              % ne×2
    bdEdge2elem = T.bdEdge2elem;    % ne×1
    clear T;
    ne = size(bdEdge,1);

    if strcmp(options.mode, 'node')
        % ---- 节点输出模式 ----
        xi = zeros(ne,3);
        bdNodeId = zeros(ne,3);
        for i = 1:ne
            [xi(i,:), bdNodeId(i,:)] = compute_on_edge_nodes(i, bdEdge, bdEdge2elem, p, t, uh);
        end
        varargout = {xi, bdEdge, bdNodeId};

    elseif strcmp(options.mode, 'gauss')
        % ---- 高斯点输出模式 ----
        nGauss = options.nGauss;
        [gp_ref, gw_ref] = gauss_rule_1d(nGauss);  % 参考区间 [-1,1]
        nTot = ne * nGauss;
        gauss_x = zeros(nTot,2);
        gauss_dudn = zeros(nTot,1);
        gauss_w = zeros(nTot,1);
        gauss_edgeIdx = zeros(nTot,1);
        idx = 0;
        for i = 1:ne
            [x_edge, dudn_gp, w_gp] = compute_on_edge_gauss(i, bdEdge, bdEdge2elem, p, t, uh, gp_ref, gw_ref);
            npt = length(dudn_gp);
            gauss_x(idx+1:idx+npt,:) = x_edge;
            gauss_dudn(idx+1:idx+npt) = dudn_gp;
            gauss_w(idx+1:idx+npt) = w_gp;
            gauss_edgeIdx(idx+1:idx+npt) = i;
            idx = idx + npt;
        end
        gauss = struct('x', gauss_x, 'dudn', gauss_dudn, 'weight', gauss_w, 'edgeIdx', gauss_edgeIdx);
        varargout = {gauss};
    end
end

% -------------------------------------------------------------------------
% 节点输出：计算边界边上三个节点的法向导数
function [xi_row, node_id_row] = compute_on_edge_nodes(edgeIdx, bdEdge, bdEdge2elem, p, t, uh)
    v1 = bdEdge(edgeIdx,1);
    v2 = bdEdge(edgeIdx,2);
    K = bdEdge2elem(edgeIdx);
    nodesK = t(K,:);
    xyK = p(nodesK,:);

    % 确定边号及边上三个节点（起点、中点、终点）
    [node_order, ~] = get_edge_nodes_order(v1, v2, nodesK);
    sample_xy = p(node_order,:);   % 3×2
    uK = uh(nodesK);
    n_outer = outer_normal(xyK(1:3,:), p([v1;v2],:));

    xi_row = zeros(1,3);
    for j = 1:3
        grad_u = compute_gradient_p2(sample_xy(j,:), xyK, uK);
        xi_row(j) = dot(grad_u, n_outer);
    end
    node_id_row = node_order;
end

% -------------------------------------------------------------------------
% 高斯点输出：计算边界边上各高斯点的法向导数及积分权重
function [x_gp, dudn_gp, w_gp] = compute_on_edge_gauss(edgeIdx, bdEdge, bdEdge2elem, p, t, uh, gp_ref, gw_ref)
    v1 = bdEdge(edgeIdx,1);
    v2 = bdEdge(edgeIdx,2);
    K = bdEdge2elem(edgeIdx);
    nodesK = t(K,:);
    xyK = p(nodesK,:);
    uK = uh(nodesK);
    n_outer = outer_normal(xyK(1:3,:), p([v1;v2],:));

    % 获取边上三个节点坐标及顺序（保证方向从v1到v2）
    [node_order, ~] = get_edge_nodes_order(v1, v2, nodesK);
    xy_edge = p(node_order,:);  % 3×2，分别为起点、中点、终点

    % 在每条边上循环高斯点
    nGauss = length(gw_ref);
    x_gp = zeros(nGauss,2);
    dudn_gp = zeros(nGauss,1);
    w_gp = zeros(nGauss,1);

    % 计算边长（用于权重缩放），此处使用端点距离近似，由于二次映射，精确Jacobian需逐点计算
    % P2边为直边，长度 = |v2-v1|
    edge_vec = xy_edge(3,:) - xy_edge(1,:);
    L = norm(edge_vec);
    if L < 1e-14
        error('边界边长度为零');
    end

    for gp = 1:nGauss
        xi = gp_ref(gp);            % 在 [-1,1]
        % 将参考坐标映射到物理边（利用三个节点的二次插值）
        N1 = 0.5*xi*(xi-1);
        N2 = 1 - xi^2;
        N3 = 0.5*xi*(xi+1);
        x_phys = N1*xy_edge(1,:) + N2*xy_edge(2,:) + N3*xy_edge(3,:);
        x_gp(gp,:) = x_phys;

        % 计算该点梯度及法向导数
        grad_u = compute_gradient_p2(x_phys, xyK, uK);
        dudn_gp(gp) = dot(grad_u, n_outer);

        % 积分权重 = 高斯权重 × (映射Jacobian)
        % 对于直边，映射为线性：x = 0.5*(v1+v2) + 0.5*xi*(v2-v1)，Jacobian = L/2
        w_gp(gp) = gw_ref(gp) * (L/2);
    end
end

% -------------------------------------------------------------------------
% 获取边界边在单元内的三个节点顺序（从v1到v2）
function [node_order, edge_loc] = get_edge_nodes_order(v1, v2, nodesK)
    vK = nodesK(1:3);
    locV1 = find(vK == v1, 1);
    locV2 = find(vK == v2, 1);
    if isempty(locV1) || isempty(locV2)
        error('边界边顶点不属于指定单元');
    end

    % 确定边号（局部编号1,2,3）及对应中点
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

    % 强制顺序：起点 v1，中点，终点 v2
    node_order = [v1, mid_node, v2];
end

% -------------------------------------------------------------------------
% 计算 P2 单元内任意物理点的梯度
function grad = compute_gradient_p2(x, xy, u)
% 通过面积坐标计算 P2 单元内任意点的梯度
A = [xy(1:3,1)'; xy(1:3,2)'; 1 1 1];
b = [x(1); x(2); 1];
L = A \ b;          % L = [L1; L2; L3]
[~, bvec, cvec] = HatGradients(xy(1:3,1), xy(1:3,2));
dNdL = zeros(6,3);
dNdL(1,1) = 4*L(1) - 1;
dNdL(2,2) = 4*L(2) - 1;
dNdL(3,3) = 4*L(3) - 1;
dNdL(4,1) = 4*L(2);  dNdL(4,2) = 4*L(1);
dNdL(5,2) = 4*L(3);  dNdL(5,3) = 4*L(2);
dNdL(6,3) = 4*L(1);  dNdL(6,1) = 4*L(3);
grad_x = (dNdL * bvec)' * u;
grad_y = (dNdL * cvec)' * u;
grad = [grad_x; grad_y];
end

% -------------------------------------------------------------------------
% 一维高斯积分规则（参考区间 [-1,1]）
function [x, w] = gauss_rule_1d(n)
    switch n
        case 1
            x = 0;
            w = 2;
        case 2
            x = [-1/sqrt(3); 1/sqrt(3)];
            w = [1; 1];
        case 3
            x = [-sqrt(3/5); 0; sqrt(3/5)];
            w = [5/9; 8/9; 5/9];
        case 4
            x = [-sqrt(3/7 + 2/7*sqrt(6/5)); -sqrt(3/7 - 2/7*sqrt(6/5));
                  sqrt(3/7 - 2/7*sqrt(6/5));  sqrt(3/7 + 2/7*sqrt(6/5))];
            w = [(18-sqrt(30))/36; (18+sqrt(30))/36; (18+sqrt(30))/36; (18-sqrt(30))/36];
        otherwise
            error('仅支持 n=1,2,3,4 的高斯规则');
    end
end