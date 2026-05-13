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