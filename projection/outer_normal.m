function n_out = outer_normal(triangle_vertices, edge_endpoints)
% outer_normal 计算三角形指定边上的单位外法向量
% 输入：
%   triangle_vertices : 3×2 矩阵，每行一个顶点坐标 (x,y)
%   edge_endpoints    : 2×2 矩阵，每行一个端点坐标 (x,y)，顺序任意
% 输出：
%   n_out : 2×1 单位外法向量

    % 提取三个顶点
    A = triangle_vertices(1,:);
    B = triangle_vertices(2,:);
    C = triangle_vertices(3,:);
    
    % 提取边的两个端点
    P1 = edge_endpoints(1,:);
    P2 = edge_endpoints(2,:);
    
    % 找出不在边上的第三个顶点
    verts = [A; B; C];
    % 判断每个顶点是否与 P1 或 P2 重合（考虑浮点误差）
    tol = 1e-12;
    idx = true(3,1);
    for i = 1:3
        if norm(verts(i,:) - P1) < tol || norm(verts(i,:) - P2) < tol
            idx(i) = false;
        end
    end
    third = verts(idx, :);  % 第三个顶点坐标
    if size(third,1) ~= 1
        error('输入的边端点不属于该三角形的边');
    end
    
    % 计算边的方向向量及单位方向
    d = P2 - P1;
    L = norm(d);
    if L < tol
        error('边长度为零');
    end
    d_unit = d / L;
    
    % 旋转90度得到法向量（逆时针旋转90度）
    n = [-d_unit(2), d_unit(1)];
    
    % 计算第三个顶点到其中一个端点（P1）的向量（顶点指向端点）
    v = P1 - third;
    
    % 点积判断方向
    if dot(n, v) > 0
        n_out = n(:);
    else
        n_out = -n(:);
    end
end