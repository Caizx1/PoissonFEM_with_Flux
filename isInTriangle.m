function z = isInTriangle(p,x,y)
% z = isInTriangle(p,x,y)
% 判断点是否在三角形中，通过计算面积实现
% 输入：
%   p：点坐标
%   x：三角形三点的横坐标
%   y：三角形三点的纵坐标
% 输出：
%   z：逻辑值，点在三角形中返回true

z = false;
p = p(:);
tol = 1e-10;

if p(1) < min(x) - tol || p(1) > max(x) + tol ||...
        p(2) < min (y) - tol || p(2) > max(y) + tol
    return;
end

area = abs(polyarea(x,y));

p1 = [x(1);y(1)] - p;
p2 = [x(2);y(2)] - p;
p3 = [x(3);y(3)] - p;
area1 = abs(det([p1,p2])) / 2;
area2 = abs(det([p2,p3])) / 2; 
area3 = abs(det([p3,p1])) / 2;
totalArea = area1 + area2 + area3;

if abs(totalArea - area) < tol
    z = true;
end