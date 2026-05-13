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