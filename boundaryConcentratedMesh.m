function [p,e,elem] = boundaryConcentratedMesh(geom,hInit,N,C,max_iter)
% 生成边界集中网格
% 输入：
%     geom：区域几何矩阵
%     hInit：初始粗网格尺寸
%     N：全局细化次数
%     C：判据常数
%     max_iter：局部细化的最大迭代次数
% 输出：
%     p：节点坐标，np * 2
%     e：边界节点信息，ne * 2
%     elem：三角形单元，nt * 3

% 初始粗网格
[p,e,elem] = initmesh(geom,"Hmax",hInit);
p = p';
elem = elem(1:3,:)';

% 保存边界线段信息，用于计算距离
bd_edges = e(1:2,:);
bd_seg = zeros(size(bd_edges,2),2,2);
for i = 1:size(bd_edges,2)
    bd_seg(i,1,:) = p(bd_edges(1,i),:);
    bd_seg(i,2,:) = p(bd_edges(2,i),:);
end

% 全局细化
for i = 1:N
    [p,elem] = bisect(p,elem);
end
h = max(elem_diam(p,elem));         % 内部网格实际尺寸

% 局部细化
for iter = 1:max_iter
   
    hT = elem_diam(p,elem);             % 单元的最大边长
    rho = dist2bd(p,elem,bd_seg);       % 单元到边界的距离

    % 标记待细化单元
    mark_elem = false(size(elem,1),1);
    for i = 1:size(mark_elem,1)
        if rho(i) < 1e-12
            if hT(i) > C * h^2
                mark_elem(i) = true;
            end
        % else
        %     if hT(i) > C * h * sqrt(rho(i))
        %         mark_elem(i) = true;
        %     end
        end

        
    end

    if ~any(mark_elem)
        break;
    end

    % 细化标记单元
    [p,elem] = bisect(p,elem,find(mark_elem));

end

% 提取边界边信息
T = auxstructure(elem);
e = T.bdEdge';

end



%------------------辅助函数------------------

function hT = elem_diam(p,elem)
% 计算单元最大边长

nT = size(elem,1);
hT = zeros(nT,1);
for i = 1:nT
    v = p(elem(i,:),:);
    l12 = norm(v(1,:)-v(2,:));
    l23 = norm(v(2,:)-v(3,:));
    l31 = norm(v(3,:)-v(1,:));
    hT(i) = max([l12,l23,l31]);
end
end

function rho = dist2bd(p,elem,bd_seg)
% 计算单元到边界的距离

nt = size(elem,1);
rho = zeros(nt,1);
tol = 1e-12;

% 检查单元是否与边界相交
for k = 1:nt
    is_intersect = false;
    nodes = elem(k,:);
    for v = 1:3
        node_coord = p(nodes(v),:);
        for j = 1:size(bd_seg,1)
            a = squeeze(bd_seg(j,1,:))';
            b = squeeze(bd_seg(j,2,:))';
            d = point2seg_dist(node_coord,a,b);
            if d < tol
                is_intersect = true;
                break;
            end
        end
        if is_intersect; break; end
    end

    if is_intersect
        rho(k) = 0;
        continue;
    end

    % 单元的重心
    centroid = (p(elem(k,1),:) + p(elem(k,2),:) + p(elem(k,3),:)) / 3;
    min_dist = inf;
    for j = 1:size(bd_seg,1)
        a = squeeze(bd_seg(j,1,:))';
        b = squeeze(bd_seg(j,2,:))';
        d = point2seg_dist(centroid,a,b);
        if d < min_dist
            min_dist = d;
        end
    end
    rho(k) = min_dist;
end
end


function d = point2seg_dist(p,a,b)
% 计算p到线段ab的距离

ab = b - a;
ap = p - a;
t = dot(ap,ab) / dot(ab,ab);

if t <= 0
    d = norm(ap);
elseif t >= 1
    d = norm(p - b);
else
    d = norm(ap - t * ab);
end
end