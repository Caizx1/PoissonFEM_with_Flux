function [B,G] = BoundaryCoupling2D_quadratic(p,t,e_b_nodes,e_boundary,g)
% [B,G] = BoundaryCoupling2D_quadratic(p,t,e_b_nodes,e_boundary,g)
% 计算原始-混合格式的边界耦合矩阵 B 和右端项 G
% 内部网格为6节点二次三角形
% 边界网格为分段线性连续元
% 输入：
%   p：内部网格节点矩阵np*2
%   t：内部网格单元nt*6
%   e_b_nodes：独立边界网格节点坐标2*nb
%   e_boundary：边界连接2*ne
%   g:边界数据函数句柄g(x,y)
% 输出：
%   B：边界耦合矩阵
%   G：右端项

np = size(p,1);
nb = size(e_b_nodes,2);
ne = size(e_boundary,2);

B = sparse(np,nb);
G = zeros(nb,1);

% 提取内部网格边界信息
T = auxstructure(t(:,1:3));
bdEdges = double(T.bdEdge);           % Nbd*2 边界边
bdEdge2elem = double(T.bdEdge2elem);
Nbd = size(bdEdges,1);

% Gauss积分参数
gp = [-1/sqrt(3),1/sqrt(3)];
gw = [1,1];

% 遍历独立边界边
tol = 1e-12;

for i = 1:ne
    idxA = e_boundary(1,i);
    idxB = e_boundary(2,i);
    pA = e_b_nodes(:,idxA);
    pB = e_b_nodes(:,idxB);
    len = norm(pB - pA);

    % 遍历内部边界边
    for k = 1:Nbd
        idxC = bdEdges(k,1);
        idxD = bdEdges(k,2);
        pC = p(idxC,:)';
        pD = p(idxD,:)';
        [t1,t2,~,~] = segment_overlap(pA,pB,pC,pD,tol);
        if t2 - t1 < tol 
            continue;
        end

        % 获取该边界所属单元
        elem_idx = bdEdge2elem(k);
        nodes = t(elem_idx,:);
        % 子区间上积分
        t_mid = (t1 + t2) / 2;
        t_half = (t2 - t1) / 2;
        for q = 1:2
            tp = t_half * gp(q) + t_mid;     % 积分点独立边上的参数
            % 物理坐标
            xq = (1 - tp) * pA(1) + tp * pB(1);
            yq = (1 - tp) * pA(2) + tp * pB(2);
            % 边界线性形函数值
            psi1 = 1 - tp;
            psi2 = tp;
            % 二次形函数
            N = quadratic_shape_functions(p,nodes,xq,yq);
            % 权重
            weight = gw(q) * t_half * len;
            % 更新B矩阵
            B(nodes,idxA) = B(nodes,idxA) + N * psi1 * weight;
            B(nodes,idxB) = B(nodes,idxB) + N * psi2 * weight;
            % 更新G矩阵
            G(idxA) = G(idxA) + weight * g(xq,yq) * psi1;
            G(idxB) = G(idxB) + weight * g(xq,yq) * psi2;
        end
    end
end

end

%-----------辅助函数------------
function N = quadratic_shape_functions(p,nodes,x,y)
% 计算二次形函数
% 输入：
%   p：网格节点坐标np*2
%   nodes：单元节点索引
%   x,y：求值点的坐标
% 输出：
%   N：形函数值向量6*1

v1 = p(nodes(1),:);
v2 = p(nodes(2),:);
v3 = p(nodes(3),:);

A = [v1(1),v2(1),v3(1);...
     v1(2),v2(2),v3(2);...
     1,    1,    1];
b = [x;y;1];
L = A \ b;

% 二次形函数
N = zeros(6, 1);
N(1) = L(1) * (2*L(1) - 1);
N(2) = L(2) * (2*L(2) - 1);
N(3) = L(3) * (2*L(3) - 1);
N(4) = 4 * L(1) * L(2);
N(5) = 4 * L(2) * L(3);
N(6) = 4 * L(3) * L(1);
end