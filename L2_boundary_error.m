function err_L2 = L2_boundary_error(e_b_nodes, e_boundary, xi, xi_exact)
% e_b_nodes: 2×nb 边界节点坐标
% e_boundary: 2×ne 边连接（每列一条边的两个端点索引）
% xi: 数值解（边界节点值，线性元）或边值（常数元）
% xi_exact: 精确函数句柄

    nb = size(e_b_nodes,2);
    ne = size(e_boundary,2);
    err = 0;

    % 高斯积分点（每条边用2点高斯，精确到线性）
    gauss_pts = [-1/sqrt(3), 1/sqrt(3)];
    gauss_w   = [1, 1];

    for i = 1:ne
        idx1 = e_boundary(1,i);
        idx2 = e_boundary(2,i);
        p1 = e_b_nodes(:,idx1);
        p2 = e_b_nodes(:,idx2);
        L = norm(p2 - p1);

        % 对每个高斯点
        for q = 1:length(gauss_pts)
            xi_g = gauss_pts(q);
            t = (xi_g+1)/2;
            x = (1-t)*p1(1) + t*p2(1);
            y = (1-t)*p1(2) + t*p2(2);
            % 数值解（线性插值）
            xi_h = (1-t)*xi(idx1) + t*xi(idx2);
            % 精确解
            xi_ex = xi_exact(x,y);
            % 误差
            err = err + (xi_h - xi_ex)^2 * (L/2) * gauss_w(q);
        end
    end
    err_L2 = sqrt(err);
end