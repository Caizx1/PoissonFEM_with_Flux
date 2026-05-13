function [B, G] = BoundaryCoupling2D_discontinuous_linear(p, t, e_b_nodes, e_boundary, g)
% [B, G] = BoundaryCoupling2D_discontinuous_linear(p, t, e_b_nodes, e_boundary, g)
% 计算原始-混合格式的边界耦合矩阵 B 和右端项 G
% 内部网格：6节点二次三角形
% 边界网格：分段线性间断元，自由度取在每个单元内部的两个高斯点
% 支持非匹配网格（边界边与内部边不重合）
%
% 输入：
%   p         : np×2 内部网格节点坐标（每行一个节点）
%   t         : nt×6 内部网格单元连接矩阵（二次三角形）
%   e_b_nodes : 2×nb 边界独立网格节点坐标（每列一个节点）
%   e_boundary: 2×ne 边界单元连接（每列一条边的两个端点索引，逆时针顺序）
%   g         : 边界数据函数句柄 g(x,y)，返回标量值
% 输出：
%   B : np × (2*ne) 稀疏矩阵，每列对应一个边界自由度的基函数
%   G : (2*ne)×1 右端项向量
%
% 边界单元上的两个高斯点（参数 t ∈ [0,1]）：
%   t1 = 0.5 - 0.5/√3， t2 = 0.5 + 0.5/√3
% 对应的间断线性基函数：
%   ψ1(t) = (t - t2)/(t1 - t2)， ψ2(t) = (t - t1)/(t2 - t1)

    % ----- 1. 提取内部网格的边界边（顶点 + 边中点）-----
    boundary_edges = extract_boundary_edges(p, t);
    n_bd_edges = length(boundary_edges);   % 内部边界边的数量

    % ----- 2. 预处理边界网格信息 -----
    np = size(p, 1);
    ne = size(e_boundary, 2);              % 边界单元数量
    nb_dof = 2 * ne;                       % 边界总自由度数
    B = sparse(np, nb_dof);
    G = zeros(nb_dof, 1);

    % 边界单元的两个高斯点（参数 t ∈ [0,1]）
    t_gauss = [0.5 - 0.5/sqrt(3), 0.5 + 0.5/sqrt(3)];
    % 参考区间 [-1,1] 上的高斯点与权重（用于子区间积分）
    tau_gauss = [-1/sqrt(3), 1/sqrt(3)];
    w_gauss = [1, 1];

    % 容差（用于 segment_overlap）
    tol = 1e-12;

    % ----- 3. 遍历每个边界单元 -----
    for e = 1:ne
        % 边界单元端点索引（注意 e_b_nodes 是 2×nb，索引从1开始）
        idxA = e_boundary(1, e);
        idxB = e_boundary(2, e);
        pA = e_b_nodes(:, idxA);
        pB = e_b_nodes(:, idxB);
        len = norm(pB - pA);

        % 边界单元对应的两列自由度的全局列号
        col1 = 2*(e-1) + 1;
        col2 = 2*(e-1) + 2;

        % 遍历所有内部边界边
        for k = 1:n_bd_edges
            % 获取内部边的节点索引和坐标
            A_int = boundary_edges(k).nodes(1);
            B_int = boundary_edges(k).nodes(2);
            M_int = boundary_edges(k).nodes(3);
            pC = p(A_int, :)';
            pD = p(B_int, :)';
            % pM 暂不需要直接使用，会在形函数计算时通过投影得到

            % 调用 segment_overlap 判断重叠区间
            [t1, t2, ~, ~] = segment_overlap(pA, pB, pC, pD, tol);
            if t2 - t1 < tol
                continue;   % 无重叠
            end

            % 子区间中点及半长（t 参数）
            t_mid = (t1 + t2) / 2;
            t_half = (t2 - t1) / 2;

            % 在重叠子区间 [t1, t2] 上做高斯积分
            for q = 1:2
                tau = tau_gauss(q);
                w = w_gauss(q);
                t = t_mid + t_half * tau;          % 边界单元局部参数 t ∈ [0,1]
                % 物理坐标
                xq = (1 - t) * pA(1) + t * pB(1);
                yq = (1 - t) * pA(2) + t * pB(2);

                % 计算内部边上三个节点的二次形函数值
                [phiA, phiB, phiM] = quadratic_shape_values_on_edge(...
                    xq, yq, pC, pD, p(A_int,:)', p(B_int,:)', p(M_int,:)');

                % 计算边界间断线性基函数值
                psi1 = (t - t_gauss(2)) / (t_gauss(1) - t_gauss(2));
                psi2 = (t - t_gauss(1)) / (t_gauss(2) - t_gauss(1));

                % 积分权重： ds = len * dt， dt = t_half * dτ，故权重 = len * t_half * w
                weight = len * t_half * w;

                % 更新 B 矩阵
                % 对于自由度 col1 (ψ1)
                B(A_int, col1) = B(A_int, col1) + phiA * psi1 * weight;
                B(B_int, col1) = B(B_int, col1) + phiB * psi1 * weight;
                B(M_int, col1) = B(M_int, col1) + phiM * psi1 * weight;
                % 对于自由度 col2 (ψ2)
                B(A_int, col2) = B(A_int, col2) + phiA * psi2 * weight;
                B(B_int, col2) = B(B_int, col2) + phiB * psi2 * weight;
                B(M_int, col2) = B(M_int, col2) + phiM * psi2 * weight;

                % 更新 G 矩阵（右端项）：∫ g ψ_k ds
                G(col1) = G(col1) + g(xq, yq) * psi1 * weight;
                G(col2) = G(col2) + g(xq, yq) * psi2 * weight;
            end
        end
    end
end

% ---------------------------------------------------------------------
% 辅助函数：提取二次三角形网格的边界边（每个边界边包含两个顶点和一个边中点）
% ---------------------------------------------------------------------
function boundary_edges = extract_boundary_edges(p, t)
    % p  : np×2 节点坐标
    % t  : nt×6 单元连接矩阵（前3列为顶点，后3列为边中点，顺序：1-2,2-3,3-1）
    % 输出 boundary_edges : 结构体数组，字段 .nodes = [A, B, M] (三个索引)

    nt = size(t, 1);
    % 收集所有边（用顶点对标识，并记录中点索引）
    % 每条边的唯一标识为排序后的两个顶点索引 [min,max]
    edges = cell(nt*3, 1);   % 预分配
    idx = 0;
    for K = 1:nt
        v1 = t(K,1); v2 = t(K,2); v3 = t(K,3);
        m12 = t(K,4); m23 = t(K,5); m31 = t(K,6);
        % 边 1-2
        idx = idx + 1;
        edges{idx} = struct('vertices', sort([v1, v2]), 'midpoint', m12, 'elem', K);
        % 边 2-3
        idx = idx + 1;
        edges{idx} = struct('vertices', sort([v2, v3]), 'midpoint', m23, 'elem', K);
        % 边 3-1
        idx = idx + 1;
        edges{idx} = struct('vertices', sort([v3, v1]), 'midpoint', m31, 'elem', K);
    end
    edges = edges(1:idx);

    % 统计每条边出现的次数
    edge_keys = cellfun(@(s) mat2str(s.vertices), edges, 'UniformOutput', false);
    [unique_keys, ~, ic] = unique(edge_keys);
    counts = accumarray(ic, 1);

    % 只保留出现一次（边界边）的边
    boundary_edges = struct('nodes', cell(0));
    for i = 1:length(unique_keys)
        if counts(i) == 1
            % 找到对应的边结构（取第一个匹配项）
            match_idx = find(strcmp(edge_keys, unique_keys{i}), 1);
            nodes = [edges{match_idx}.vertices(1), edges{match_idx}.vertices(2), edges{match_idx}.midpoint];
            boundary_edges(end+1).nodes = nodes;
        end
    end
end

% ---------------------------------------------------------------------
% 辅助函数：计算边上一点处三个二次形函数的值（顶点A，顶点B，边中点M）
% 基于参数 t ∈ [0,1]，其中 t=0 对应 A，t=1 对应 B
% ---------------------------------------------------------------------
function [phiA, phiB, phiM] = quadratic_shape_values_on_edge(x, y, pA, pB, pA_coord, pB_coord, pM_coord)
    % 输入：
    %   x,y       : 物理坐标
    %   pA, pB    : 端点坐标（列向量，2×1）
    %   pA_coord, pB_coord, pM_coord : 三个节点坐标（行向量，1×2）
    % 实际上 pA 和 pA_coord 相同，此处为了接口统一传递即可。
    % 计算参数 t (投影到 AB 上)
    AB = pB - pA;
    len2 = AB' * AB;
    if len2 < eps
        t = 0;
    else
        t = ((x - pA(1)) * AB(1) + (y - pA(2)) * AB(2)) / len2;
        t = max(0, min(1, t));   % 确保 t 在 [0,1] 内（浮点误差）
    end
    % 二次形函数公式（参考区间 [0,1]）
    phiA = (1 - t) * (1 - 2*t);
    phiB = t * (2*t - 1);
    phiM = 4 * t * (1 - t);
end