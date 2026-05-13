function err = H1_error(p, t, u_h, u_exact, grad_u_exact, order)
    % err = H1_error(p,t,u_h,u_exact,grad_u_exact,order)
    % 计算H1范数下的误差
    % p: np×2 节点坐标
    % t: nt×3 (linear) 或 nt×6 (quadratic)
    % u_h: 节点解 (np×1)
    % u_exact: 精确解函数句柄
    % grad_u_exact: 精确梯度函数句柄，返回 [∂u/∂x, ∂u/∂y]
    % order: 'linear' 或 'quadratic'
    
    nt = size(t,1);
    err = 0;
    
    % 高斯积分点（面积坐标）和权重（三点积分，精确到5次）
    gp = [1/6, 1/6, 2/3, 1/6;
          1/6, 2/3, 1/6, 1/6;
          2/3, 1/6, 1/6, 1/6];
    nGauss = size(gp,1);
    
    for K = 1:nt
        if strcmp(order, 'linear')
            nodes = t(K,1:3);
            nNodes = 3;
        else
            nodes = t(K,:);
            nNodes = 6;
        end
        xy = p(nodes, :);   % nNodes × 2
        ue = u_h(nodes);    % nNodes × 1
        
        errK = 0;
        for i = 1:nGauss
            L1 = gp(i,1); L2 = gp(i,2); L3 = gp(i,3); w = gp(i,4);
            if strcmp(order, 'linear')
                % 线性形函数（面积坐标）
                N = [L1, L2, L3];   % 1×3
                % 形函数对 L1, L2 的导数（对 L3 = 1-L1-L2）
                dNdL = [1, 0; 0, 1; -1, -1];  % 3×2
            else
                % 二次形函数
                N = [L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1), 4*L1*L2, 4*L2*L3, 4*L3*L1]; % 1×6
                % 对 L1, L2 的偏导数（6×2）
                dNdL = zeros(6,2);
                dNdL(1,1) = 4*L1 - 1;  dNdL(1,2) = 0;
                dNdL(2,1) = 0;         dNdL(2,2) = 4*L2 - 1;
                dNdL(3,1) = 1 - 4*L3;  dNdL(3,2) = 1 - 4*L3;
                dNdL(4,1) = 4*L2;      dNdL(4,2) = 4*L1;
                dNdL(5,1) = -4*L2;     dNdL(5,2) = 4*L3 - 4*L2;
                dNdL(6,1) = 4*L3 - 4*L1; dNdL(6,2) = -4*L1;
            end
            % Jacobian 矩阵 (2×2)
            J = xy' * dNdL;   % xy' 是 2×nNodes，dNdL 是 nNodes×2，乘积 2×2
            detJ = abs(det(J));
            % 物理坐标
            xq = N * xy(:,1);   % 点积：1×nNodes * nNodes×1 = 标量
            yq = N * xy(:,2);
            % 形函数对物理坐标的导数 (nNodes×2)
            dNdx = dNdL / J;   
            % 数值解和梯度
            u_h_q = N * ue;          % 1×nNodes * nNodes×1 = 标量
            grad_u_h_q = ue' * dNdx;  % 1×nNodes * nNodes×2 = 1×2
            % 精确解和梯度
            u_exact_q = u_exact(xq, yq);
            grad_u_exact_q = grad_u_exact(xq, yq);
            % 误差
            e = u_exact_q - u_h_q;
            grad_e = grad_u_exact_q - grad_u_h_q;
            % 累加
            errK = errK + (e^2 + grad_e(1)^2 + grad_e(2)^2) * w * detJ;
        end
        err = err + errK;
    end
    err = sqrt(err);
end