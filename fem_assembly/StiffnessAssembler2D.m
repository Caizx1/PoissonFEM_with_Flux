function A = StiffnessAssembler2D(p,t,a,order)
% 组装刚度矩阵
% 输入：
%   p：2 * np，节点坐标
%   t：4 * nt（线性单元）或7 * nt（二次单元），单元连接矩阵
%   a：系数函数，函数句柄@(x,y)
%   order：基函数阶数，可选'linear'（默认）或'quadratic'
% 输出：
%   A：整体刚度矩阵

if nargin < 4
    order = 'linear';
end

np = size(p,2);
nt = size(t,2);
A = sparse(np,np);

if strcmp(order,'linear')
    % 线性单元
    for K = 1:nt
        loc2glb = t(1:3,K);
        x = p(1,loc2glb);
        y = p(2,loc2glb);
        [area,b,c] = HatGradients(x,y);
        xc = mean(x);yc = mean(y);      %单元重心
        abar = a(xc,yc);                %单元刚性矩阵，重心积分公式
        AK = abar * (b * b' + c * c') * area;
        A(loc2glb,loc2glb) = A(loc2glb,loc2glb) + AK;
    end

elseif strcmp(order,'quadratic')
    % 二次单元
    % 三点高斯积分（面积坐标）
    % 积分点: (L1,L2,L3) 和权重 w
    gp = [1/6, 1/6, 2/3, 1/6;
          1/6, 2/3, 1/6, 1/6;
          2/3, 1/6, 1/6, 1/6];
    nGauss = size(gp,1);

    for K = 1:nt
        loc2glb = t(1:6,K);
        xy = p(:,loc2glb);
        Ke = zeros(6,6);        % 单元刚度矩阵

        for i = 1:nGauss
            L1 = gp(i,1); L2 = gp(i,2); L3 = gp(i,3); w = gp(i,4);
            % 形函数（1*6）
            N = [L1 * (2 * L1 - 1),L2 * (2 * L2 - 1),L3 * (2 * L3 - 1),...
                4 * L1 * L2,4 * L2 * L3,4 * L3 * L1];
            % 形函数对L1和L2的偏导数（6*2）
            dNdL = [4 * L1 - 1,     0;...
                    0,              4 * L2 - 1;...
                    1 - 4 * L3,     1 - 4 * L3;...
                    4 * L2,         4 * L1;...
                    -4 * L2,        4 * L3 - 4 * L2;...
                    4 * L3 - 4 * L1,-4 * L1];
            % Jacobian矩阵（2*2）
            J = xy * dNdL;
            detJ = abs(det(J));
            % 形函数对物理坐标的导数（6*2）
            dNdx = dNdL / J;
            % 高斯点物理坐标
            xp = xy(1,:) * N';
            yp = xy(2,:) * N';

            Ke = Ke + w * a(xp,yp) * detJ * (dNdx * dNdx') ;
        end
        A(loc2glb,loc2glb) = A(loc2glb,loc2glb) + Ke;
    end

else
    error('order参数须是"linear"或"quadratic"');
end

end