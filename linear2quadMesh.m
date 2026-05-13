function [p_quad,t_quad] = linear2quadMesh(p_lin,t_lin)
% 将线性三角形网格转换为二次网格（6节点）
% 输入：
%   p_lin：节点坐标np*2
%   t_lin：连接矩阵nt*3
% 输出：
%   p_quad：节点坐标(np + n_edges)*2
%   t_quad：连接矩阵nt*6

np = size(p_lin,1);
nt = size(t_lin,1);

% 收集所有边
edges = [t_lin(:,[1 2]);t_lin(:,[2 3]);t_lin(:,[3 1])];
edges = sort(edges,2);
[edges,~,idx_e] = unique(edges,"rows");

% 单元的边在edges中的序号
elem2edge = reshape(idx_e,nt,3);

% 更新节点
midpoint = (p_lin(edges(:,1),:) + p_lin(edges(:,2),:)) / 2;
p_quad = [p_lin;midpoint];

% 更新连接矩阵
t_quad = zeros(nt,6);
for K = 1:nt
    mid12 = elem2edge(K,1) + np;
    mid23 = elem2edge(K,2) + np;
    mid31 = elem2edge(K,3) + np;

    t_quad(K,:) = [t_lin(K,:),mid12,mid23,mid31];
end