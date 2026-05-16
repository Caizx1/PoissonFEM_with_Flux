function [pnew, tnew] = trimeshRefine(p, t)
% 均匀红绿细化：将每个三角形一分为四（通过边中点连接）
% 输入：
%   p - np×2 节点坐标
%   t - nt×3 三角形单元（每行三个顶点索引）
% 输出：
%   pnew - (np + n_edges)×2 新节点坐标
%   tnew - (4*nt)×3 新三角形单元

    np = size(p, 1);
    nt = size(t, 1);

    % 1. 收集所有边（单元的三条边）
    edges = [t(:, [1,2]); t(:, [2,3]); t(:, [3,1])];   % 3*nt × 2
    edges = sort(edges, 2);   % 每行从小到大
    [edges,~,idx_e] = unique(edges, 'rows');   % 去重，每行一条边（顶点1,顶点2）

    % 2. 计算每条边的中点，并添加到节点列表
    midpoints = (p(edges(:,1), :) + p(edges(:,2), :)) / 2;  % n_edges × 2
    pnew = [p; midpoints];    % 新节点列表，前 np 个为原节点，后 n_edges 个为中点

    % 3. 建立单元与边的映射
    elem2edge = reshape(idx_e,nt,3);

    % 4. 细化每个三角形
    tnew = zeros(4*nt, 3);
    for K = 1:nt
        v1 = t(K,1); v2 = t(K,2); v3 = t(K,3);

        m12 = elem2edge(K,1) + np;
        m23 = elem2edge(K,2) + np;
        m31 = elem2edge(K,3) + np;

        idx = (K-1)*4 + 1;
        tnew(idx,   :) = [v1, m12, m31];
        tnew(idx+1, :) = [v2, m23, m12];
        tnew(idx+2, :) = [v3, m31, m23];
        tnew(idx+3, :) = [m12, m23, m31];
    end
end