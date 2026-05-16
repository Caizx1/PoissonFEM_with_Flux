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