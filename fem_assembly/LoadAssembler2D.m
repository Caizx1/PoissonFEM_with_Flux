function b = LoadAssembler2D(p,t,f,order)
% 组装载荷矩阵
% 输入：
%   p：2*np，节点坐标
%   t：4*nt（线性）或7*nt（二次），单元连接矩阵
%   f：源项函数，函数句柄f(x,y)
%   order：基函数阶数，可选'linear'（默认）或'quadratic'
% 输出
%   b：载荷矩阵

if nargin < 4
    order = 'linear';
end

np = size(p,2);
nt = size(t,2);
b = sparse(np,1);

if strcmp(order,'linear')
    % 线性单元，使用节点平均值近似积分
    for K = 1:nt
        loc2glb = t(1:3,K);
        x = p(1,loc2glb);
        y = p(2,loc2glb);
        area = polyarea(x,y);
        bK = [f(x(1),y(1));f(x(2),y(2));f(x(3),y(3))] / 3 * area;
        b(loc2glb) = b(loc2glb) + bK;
    end

elseif strcmp(order,'quadratic')
    % 二次单元，三点高斯积分
    % 积分点: (L1,L2,L3) 和权重 w
    gp = [1/6, 1/6, 2/3, 1/6;
          1/6, 2/3, 1/6, 1/6;
          2/3, 1/6, 1/6, 1/6];
    nGauss = size(gp,1);

    for K = 1:nt
        loc2glb = t(1:6,K);
        xy = p(:,loc2glb);
        be = zeros(6,1);

        for i = 1:nGauss
            L1 = gp(i,1); L2 = gp(i,2); L3 = gp(i,3); w = gp(i,4);
            % 形函数（1*6）
            N = [L1 * (2 * L1 - 1),L2 * (2 * L2 - 1),L3 * (2 * L3 - 1),...
                4 * L1 * L2,4 * L2 * L3,4 * L3 * L1];
            % 形函数对L1和L2的偏导数（2*6）
            dNdL = [4 * L1 - 1,     0;...
                    0,              4 * L2 - 1;...
                    1 - 4 * L3,     1 - 4 * L3;...
                    4 * L2,         4 * L1;...
                    -4 * L2,        4 * L3 - 4 * L2;...
                    4 * L3 - 4 * L1,-4 * L1]';
            % Jacobian矩阵（2*2）
            J = xy * dNdL';
            detJ = abs(det(J));
            % 高斯点物理坐标
            xp = xy(1,:) * N';
            yp = xy(2,:) * N';

            be = be + w * detJ * f(xp,yp) * N';
        end
        b(loc2glb) = b(loc2glb) + be;
    end

else
    error('order须是"linear"或"quadratic"');
end
end