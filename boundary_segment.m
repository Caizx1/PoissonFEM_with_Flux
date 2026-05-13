function bd_seg = boundary_segment(geom)
% 从区域几何矩阵中提取边界信息

n = size(geom,2);
bd_seg = zeros(n,2,2);
for i = 1:n
    bd_seg(i,1,:) = geom([2,4],i)';
    bd_seg(i,2,:) = geom([3,5],i)';
end
end