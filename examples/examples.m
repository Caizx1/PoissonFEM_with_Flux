function examples = examples()
% examples = examples() 返回预定义的算例结构体数组
% 每个结构体包含字段：
%   name      - 算例名称
%   geom      - 几何矩阵
%   u_exact   - 精确解函数句柄
%   grad_u_exact - 梯度函数句柄
%   f         - 右端项函数句柄
%   g         - Dirichlet 边界条件函数句柄
%   xi_exact_fun - 边界法向导数精确值函数句柄（ξ = -∂u/∂n）

    % 算例1：齐次 Dirichlet，光滑解 u = sin(πx)sin(πy)
    geom = Rectg(0, 0, 1, 1);
    examples(1).name = '齐次 Dirichlet (光滑解)';
    examples(1).geom = geom;
    examples(1).u_exact = @(x,y) sin(pi*x) .* sin(pi*y);
    examples(1).grad_u_exact = @(x,y) [pi*cos(pi*x).*sin(pi*y), pi*sin(pi*x).*cos(pi*y)];
    examples(1).f = @(x,y) 2*pi^2 * sin(pi*x) .* sin(pi*y);
    examples(1).g = @(x,y) 0;
    examples(1).xi_exact_fun = @(x,y) ...
        (abs(x)<1e-12) * (-pi * sin(pi*y)) + ...
        (abs(x-1)<1e-12) * (-pi * sin(pi*y)) + ...
        (abs(y)<1e-12) * (-pi * sin(pi*x)) + ...
        (abs(y-1)<1e-12) * (-pi * sin(pi*x));

    % 算例2：非齐次 Dirichlet，解 u = exp(x)sin(πy) + exp(y)sin(πx)，角点奇异性
    examples(2).name = '非齐次 Dirichlet (指数-三角, 角点奇异性)';
    examples(2).geom = geom; % 复用几何
    examples(2).u_exact = @(x,y) exp(x).*sin(pi*y) + exp(y).*sin(pi*x);
    examples(2).grad_u_exact = @(x,y) [exp(x).*sin(pi*y) + pi*exp(y).*cos(pi*x), ...
                                      pi*exp(x).*cos(pi*y) + exp(y).*sin(pi*x)];
    examples(2).f = @(x,y) (pi^2 - 1) * (exp(x).*sin(pi*y) + exp(y).*sin(pi*x));
    examples(2).g = examples(2).u_exact;
    examples(2).xi_exact_fun = @(x,y) ...
        (abs(y)   < 1e-12) * (-( pi*exp(x).*cos(pi*y) + exp(y).*sin(pi*x) )) + ...
        (abs(y-1) < 1e-12) * (  pi*exp(x).*cos(pi*y) + exp(y).*sin(pi*x) ) + ...
        (abs(x)   < 1e-12) * ( -(exp(x).*sin(pi*y) + pi*exp(y).*cos(pi*x)) ) + ...
        (abs(x-1) < 1e-12) * (  exp(x).*sin(pi*y) + pi*exp(y).*cos(pi*x) );

    % 算例3：非齐次 Dirichlet，解 u = sin(πx)+sin(πy)，法向导数连续
    examples(3).name = '非齐次 Dirichlet (正弦和, 法向导数连续)';
    examples(3).geom = geom;
    examples(3).u_exact = @(x,y) sin(pi*x) + sin(pi*y);
    examples(3).grad_u_exact = @(x,y) [pi*cos(pi*x), pi*cos(pi*y)];
    examples(3).f = @(x,y) pi^2 * (sin(pi*x) + sin(pi*y));
    examples(3).g = examples(3).u_exact;
    examples(3).xi_exact_fun = @(x,y) -pi * ones(size(x));


    % 算例4：精确解 u = x^2*y^2*(2*x-3)*(2*y-3) + sin(pi*x)*sin(pi*y)
    examples(4).name = '非齐次 Dirichlet (多项式+三角, 通量连续)';
    examples(4).geom = Rectg(0, 0, 1, 1);
    examples(4).u_exact = @(x,y) x.^2.*y.^2.*(2*x-3).*(2*y-3) + sin(pi*x).*sin(pi*y);
    examples(4).grad_u_exact = @(x,y) [ ...
        (6*x.*(x-1)).*(2*y.^3-3*y.^2) + pi*cos(pi*x).*sin(pi*y), ...
        (2*x.^3-3*x.^2).*(6*y.*(y-1)) + pi*sin(pi*x).*cos(pi*y) ];
    examples(4).f = @(x,y) ...
        -24*x.*y.^3 + 36*x.*y.^2 + 12*y.^3 - 18*y.^2 ...
        -24*x.^3.*y + 12*x.^3 + 36*x.^2.*y - 18*x.^2 ...
        + 2*pi^2 * sin(pi*x).*sin(pi*y);
    examples(4).g = examples(4).u_exact;   % 非齐次 Dirichlet
    examples(4).xi_exact_fun = @(x,y) ...
        (abs(x)<1e-12) .* (-pi * sin(pi*y)) + ...      % 左边界
        (abs(x-1)<1e-12) .* (-pi * sin(pi*y)) + ...    % 右边界
        (abs(y)<1e-12) .* (-pi * sin(pi*x)) + ...      % 下边界
        (abs(y-1)<1e-12) .* (-pi * sin(pi*x));         % 上边界

    % 算例5：u = x^2 + y^2 + sin(πx)sin(πy)，多项式部分在边界上法向导数不消失
    % 导致角点(0,1)和(1,0)处通量不连续，其他两角连续
    examples(5).name = '非齐次 Dirichlet (多项式+三角, 角点通量不连续)';
    examples(5).geom = Rectg(0, 0, 1, 1);
    examples(5).u_exact = @(x,y) x.^2 + y.^2 + sin(pi*x).*sin(pi*y);
    examples(5).grad_u_exact = @(x,y) [ ...
        2*x + pi*cos(pi*x).*sin(pi*y), ...
        2*y + pi*sin(pi*x).*cos(pi*y) ];
    examples(5).f = @(x,y) -4 + 2*pi^2 * sin(pi*x).*sin(pi*y);
    examples(5).g = examples(5).u_exact;
    examples(5).xi_exact_fun = @(x,y) ...
        (abs(x)<1e-12) .* (-pi * sin(pi*y)) + ...      % 左边界: ξ = -∂u/∂x = -π sin(πy)
        (abs(x-1)<1e-12) .* (2 - pi * sin(pi*y)) + ... % 右边界: ξ = ∂u/∂x = 2 - π sin(πy)
        (abs(y)<1e-12) .* (-pi * sin(pi*x)) + ...      % 下边界: ξ = -∂u/∂y = -π sin(πx)
        (abs(y-1)<1e-12) .* (2 - pi * sin(pi*x));      % 上边界: ξ = ∂u/∂y = 2 - π sin(πx)

    % 算例6：u = x + 2y + sin(πx)sin(πy)，利用不对称系数(1≠2)
    % 使四个角点的法向导数均不连续
    examples(6).name = '非齐次 Dirichlet (不对称线性项+三角, 四角点通量不连续)';
    examples(6).geom = Rectg(0, 0, 1, 1);
    examples(6).u_exact = @(x,y) x + 2*y + sin(pi*x).*sin(pi*y);
    examples(6).grad_u_exact = @(x,y) [ ...
        1 + pi*cos(pi*x).*sin(pi*y), ...
        2 + pi*sin(pi*x).*cos(pi*y) ];
    examples(6).f = @(x,y) 2*pi^2 * sin(pi*x).*sin(pi*y);
    examples(6).g = examples(6).u_exact;
    examples(6).xi_exact_fun = @(x,y) ...
        (abs(x)<1e-12) .* (-1 - pi*sin(pi*y)) + ...     % 左边界: ξ = -∂u/∂x = -1 - π sin(πy)
        (abs(x-1)<1e-12) .* (1 - pi*sin(pi*y)) + ...    % 右边界: ξ = ∂u/∂x = 1 - π sin(πy)
        (abs(y)<1e-12) .* (-2 - pi*sin(pi*x)) + ...     % 下边界: ξ = -∂u/∂y = -2 - π sin(πx)
        (abs(y-1)<1e-12) .* (2 - pi*sin(pi*x));         % 上边界: ξ = ∂u/∂y = 2 - π sin(πx)

end