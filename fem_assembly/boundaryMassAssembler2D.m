function M = boundaryMassAssembler2D(p,t,bdNodes,order)
% 计算边界上的质量矩阵
% 输入：
%   p：节点坐标矩阵2*np
%   t：连接矩阵3*nt(linear)或6*nt
%   bdNodes:边界节点向量
%   order：基函数阶数，可选'liear'（默认）和'quadratic'
% 输出：
%   M：边界质量矩阵

if nargin < 4
    order = 'linear';
end

T = auxstructure(t(1:3,:)');
bdEdgel = T.bdEdge;
n_bdEdgel = size(bdEdgel,1);
n_bdNodes = size(bdNodes,1);
M = sparse(n_bdNodes,n_bdNodes);

if strcmp(order,'linear')
    for e = 1:n_bdNodes
        % 端点索引
        n1 = bdEdgel(e,1);
        n2 = bdEdgel(e,2);
        % 端点在质量矩阵的索引
        idx1 = find(bdNodes == n1);
        idx2 = find(bdNodes == n2);
        % 端点坐标
        p1 = p(:,n1);
        p2 = p(:,n2);
        % 边长
        L = norm(p2 - p1);
    
        Me = [2,1;1,2] / 6 * L;
        M([idx1,idx2],[idx1,idx2]) = M([idx1,idx2],[idx1,idx2]) + Me;
            
    end

elseif strcmp(order,'quadratic')
    for e = 1: n_bdEdgel
        % 端点索引
        n1 = bdEdgel(e,1);
        n2 = bdEdgel(e,2);
        % 节点坐标
        p1 = p(:,n1);
        p2 = p(:,n2);
        
        pm = (p1 + p2) / 2;
        mid_idx = all(abs(p - pm) < 1e-10,1);
        n3 = find(mid_idx);         % 中点索引
        % 端点在质量矩阵中索引
        idx1 = find(bdNodes == n1);
        idx2 = find(bdNodes == n2);
        idx3 = find(bdNodes == n3);
        loc2glb = [idx1,idx2,idx3];
        
        % 边长
        L = norm(p2 - p1);

        Me = [4, -1, 2; -1, 4, 2; 2, 2, 16] / 30 * L;
        M(loc2glb,loc2glb) = M(loc2glb,loc2glb) + Me;
    end
end
        
end