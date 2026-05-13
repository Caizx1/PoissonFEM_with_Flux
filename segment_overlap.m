function [t1,t2,s1,s2] = segment_overlap(A,B,C,D,tol)
% [t1,t2,s1,s2] = segment_overlap(A,B,C,D,tol)
% 计算两条线段的重叠参数
% 输入：
%   A，B：外部边的端点坐标
%   C，D：内部边的端点坐标
%   tol:容差（默认1e-12)
% 输出：
%   t1,t2：重叠区间在 AB 上的参数（0~1）
%   s1,s2：在 CD 上的参数
%   若没有重叠区间，返回 t1=t2=0, s1=s2=0

if nargin < 5
    tol = 1e-12;
end

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