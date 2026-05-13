function [B,G] = BoundaryCoupling2D0(p,e_b,g)
%计算原始-混合格式的边界耦合矩阵

np = size(p,2);
ne_b = size(e_b,2);
B = sparse(np,ne_b);
G = sparse(ne_b,1);

for Eb = 1:ne_b
    p1 = e_b(1,Eb);
    p2 = e_b(2,Eb);
    length_eb = norm(p(:,p1) - p(:,p2));    %边的长度
    
    B(p1,Eb) = B(p1,Eb) + length_eb / 2;    
    B(p2,Eb) = B(p2,Eb) + length_eb / 2;

    G(Eb) = g((p(:,p1) + p(:,p2)) / 2) * length_eb;   %中点公式
end