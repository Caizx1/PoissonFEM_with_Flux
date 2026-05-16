function points = boundaryDivide(geom,step)
% 将区域的边界进行分割
% 输入：
%   geom：区域的几何矩阵7 * N 
%   step：最大尺寸
% 输出：
%   points：分割点的坐标矩阵2*Ne
  

N = size(geom,2);       %边界数量
points = [];            

for i = 1:N
    %端点坐标
    p1 = geom([2,4],i);
    p2 = geom([3,5],i);

    %计算线段长度
    l = norm(p2 - p1);

    n = ceil(l / step);     %划分数
    t = linspace(0,1,n + 1);
    t = t(1:n);             %取前n个点，不包含终点
    seg = p1 + t .* (p2 - p1);
    points = [points,seg];
end