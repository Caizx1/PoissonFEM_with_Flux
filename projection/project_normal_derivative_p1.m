function [xi,bdEdge] = project_normal_derivative_p1(uh,p,t)
% xi = project_normal_derivative(uh,p,e,t)
% 计算投影法向导数
% 输入：
%   uh：节点解np*1
%   p：节点坐标np*2
%   t：三角形单元nt*3
% 输出：
%   xi：边界法向导数
%   bdEdge：边界边信息ne*2

T = auxstructure(t);
bdEdge = T.bdEdge;              % 边界边ne*2
bdEdge2elem = T.bdEdge2elem;    % 边界边所在单元
clear T;

ne = size(bdEdge2elem,1);
xi = zeros(ne,1);

for i = 1:ne
    edge_end = bdEdge(i,:);
    K = bdEdge2elem(i);
    nodes = t(K,:);
    edge_end_points = p(edge_end,:);
    xy = p(nodes,:);

    [~,b,c] = HatGradients(xy(:,1),xy(:,2));
    grad_u = [b,c]' * uh(nodes);        % u的梯度2*1
    n_outer = outer_normal(xy,edge_end_points); % 单位外法向量
    dudn = grad_u' * n_outer;
    xi(i) = dudn;
end

end