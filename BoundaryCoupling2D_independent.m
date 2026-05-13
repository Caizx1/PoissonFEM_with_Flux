function [B, G] = BoundaryCoupling2D_independent(p, e, e_b, g)
% [B, G] = BoundaryCoupling2D_independent(p, e, e_b, g)
% 计算原始-混合格式的边界耦合矩阵 B 和右端项 G，支持内部网格与边界网格独立
% 输入参数：
%   p   - 内部网格节点坐标，2×np
%   e   - 内部网格边矩阵，由 initmesh 返回，格式为 7×ne
%   e_b - 边界网格节点坐标，大小 2×ne_node，按逆时针顺序排列
%   g   - 边界数据函数句柄，g(x,y) 返回标量值
% 输出参数：
%   B   - 耦合矩阵，大小 np×ne_edge
%   G   - 边界数据向量，大小 ne_edge×1

np = size(p,2);
ne_int = size(e,2);
ne_edge = size(e_b,2);
B = sparse(np,ne_edge);
G = zeros(ne_edge,1);

tol = 1e-12;               %容差

for i = 1:ne_edge
    if i < ne_edge
        i1 = i;i2 = i + 1;
    else
        i1 = ne_edge;i2 = 1;
    end

    pA = e_b(:,i1);pB = e_b(:,i2);      % 外部边端点坐标

    for k = 1:ne_int
        idxC = e(1,k);idxD = e(2,k);            % 内部边端点索引
        pC = p(:,idxC); pD = p(:,idxD);         % 内部边端点坐标
        [t1,t2,s1,s2] = segment_overlap(pA,pB,pC,pD,tol);   %判断内外边是否重叠
        if t2 - t1 < tol
            continue
        end
        len_overlap = norm(pC - pD) * (s2 - s1);            %重叠部分长度
        
        sum_s = s1 + s2;
        %更新B矩阵
        B(idxC,i) = B(idxC,i) + len_overlap * (1 - sum_s / 2);
        B(idxD,i) = B(idxD,i) + len_overlap * sum_s / 2;

        % 更新G矩阵
        mid_s = (s1 + s2) / 2;
        xm = (1 - mid_s) * pC(1) + mid_s * pD(1);
        ym = (1 - mid_s) * pC(2) + mid_s * pD(2);
        G(i) = G(i) + g(xm, ym) * len_overlap;
    end
end
    
end

%-------------------------------

function [t1,t2,s1,s2] = segment_overlap(A,B,C,D,tol)
% [t1,t2,s1,s2] = segment_overlap(A,B,C,D,tol)
% 计算两条线段的重叠参数
% 输入：
%   A，B：外部边的端点坐标
%   C，D：内部边的端点坐标
%   tol:容差
% 输出：
%   t1,t2：重叠区间在 AB 上的参数（0~1）
%   s1,s2：在 CD 上的参数
%   若没有重叠区间，返回 t1=t2=0, s1=s2=0

t1 = 0;t2 = 0;s1 = 0;s2 = 0;
A = A(:);B = B(:);C = C(:);D = D(:);

v = B - A;
w = D - C;
M = [v,-w];

% 交叉
if abs(det(M)) > tol
    return
end

% 平行但不共线
cross1 = det([v,C - A]);
if abs(cross1) > tol
    return
end

% 计算C，D相对v的参数
len2 = v' * v;
tC = (C - A)' * v / len2;
tD = (D - A)' * v / len2;

% 重叠部分在v的参数
t_start = max(0,min(tC,tD));
t_end = min(1,max(tC,tD));
% 共线无重叠
if t_end - t_start < tol
    return
end

% 重叠部分在CD上的参数
s_start = (t_start - tC) / (tD - tC);
s_end   = (t_end   - tC) / (tD - tC);
if s_start > s_end
    [s_start, s_end] = deal(s_end, s_start);
end

t1 = t_start; t2 = t_end;
s1 = s_start; s2 = s_end;

end