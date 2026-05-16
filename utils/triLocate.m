function tri = triLocate(p,t,coord,candidateElems)
% idx = triLocat(p,t,coord,candidateElems)
% 返回点所在的第一个候选单元索引
% 输入：
%   p：np*2，节点坐标
%   t：nt*3，三角单元连接矩阵
%   coord：点坐标
%   candidateElems：候选单元索引列表（默认1:nt)
% 输出：
%   tri:点所在的第一个单元的索引

if nargin < 4
    candidateElems = 1:size(t,1);
end

tri = nan;
for k = 1:length(candidateElems)
    idx = candidateElems(k);
    nodes = t(idx,:);
    xv = p(nodes,1);
    yv = p(nodes,2);
    if isInTriangle(coord,xv,yv)
        tri = idx; 
        return;
    end
end
