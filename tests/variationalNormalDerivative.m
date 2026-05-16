%使用变分法计算法向导数

clear;clc;
% ===================== 定义问题和精确解 =====================
geom = Rectg(0, 0, 1, 1);   % 几何

u_exact = @(x,y) sin(pi * x) .* sin(pi * y);
a = @(x,y)1;
f = @(x,y) 2 * pi ^ 2 * sin(pi * x) .* sin(pi * y);
g = @(x,y) 0;
xi_exact_fun = @(x,y) ...
    (abs(x)<1e-12) * (-pi) .* sin(pi * y) + ...
    (abs(x-1)<1e-12) * (-pi) .* sin(pi * y) + ...
    (abs(y)<1e-12) * (-pi) .* sin(pi * x) + ...
    (abs(y-1)<1e-12) * (-pi) .* sin(pi * x);

h = [0.2 0.1 0.05 0.025 0.0125];
err_l_L2 = [];
err_q_L2 = [];

% fprintf('使用加权L2范数计算误差\n');
fprintf('使用离散l2范数计算误差\n');
for i = 1:size(h,2)
    % 定义网格
    [pl,el,tl] = initmesh(geom,"Hmax",h(i));
    [pq,tq] = linear2quadMesh(pl',tl(1:3,:)');
    pq = pq';tq = tq';
    
    % 矩阵组装
    Al_origin = StiffnessAssembler2D(pl,tl,a);
    bl_origin = LoadAssembler2D(pl,tl,f);
    Al = Al_origin;
    bl = bl_origin;

    Aq_origin = StiffnessAssembler2D(pq,tq,a,'quadratic');
    bq_origin = LoadAssembler2D(pq,tq,f);
    Aq = Aq_origin;
    bq = bq_origin;
    
    % 提取边界信息
    Tl = auxstructure(tl(1:3,:)');
    bdEdge_l = Tl.bdEdge;
    bdNodes_l = bdEdge_l(:,1);
    bdpoints_l = pl(:,bdNodes_l);
    n_bdNodes_l = length(bdNodes_l);

    [pqn,tqn] = trimeshRefine(pl',tl(1:3,:)');
    pqn = pqn';tqn = tqn';
    Tq = auxstructure(tqn');
    bdEdge_q = Tq.bdEdge;
    bdNodes_q = bdEdge_q(:,1);
    bdpoints_q = pqn(:,bdNodes_q);
    nbdNodes_q = size(bdNodes_q,1);
    
    % 加载边界条件
    Al(:,bdNodes_l) = 0;Al(bdNodes_l,:) = 0;
    Al(bdNodes_l,bdNodes_l) = speye(n_bdNodes_l);
    bl(bdNodes_l) = g(bdpoints_l(1,:),bdpoints_l(2,:));

    Aq(bdNodes_q,:) = 0; Aq(:,bdNodes_q) = 0;
    Aq(bdNodes_q,bdNodes_q) = speye(nbdNodes_q);
    bq(bdNodes_q) = g(bdpoints_q(1,:),bdpoints_q(2,:));
    
    % 求解，计算残差
    ul = Al \ bl;
    rl = Al_origin * ul - bl_origin;
    
    Mb_l = boundaryMassAssembler2D(pl,tl,bdNodes_l);
    rb_l = rl(bdNodes_l);

    uq = Aq \ bq;
    rq = Aq_origin * uq - bq_origin;

    Mb_q = boundaryMassAssembler2D(pq,tq,bdNodes_q,'quadratic');
    rb_q = rq(bdNodes_q);
    
    % 法向导数计算
    dudn_l = Mb_l \ rb_l;
    dudn_l_exact = xi_exact_fun(pl(1,bdNodes_l),pl(2,bdNodes_l));
    dudn_l_exact = dudn_l_exact(:);

    dudn_q = Mb_q \ rb_q;
    dudn_q_exact = xi_exact_fun(pq(1,bdNodes_q),pq(2,bdNodes_q));
    dudn_q_exact = dudn_q_exact(:);

    % 误差计算
    err_l = dudn_l_exact - dudn_l;
    err_l_L2 = [err_l_L2,sqrt(err_l' * err_l / length(err_l))];
    % err_l_L2 = [err_l_L2,sqrt(err_l' * Mb_l * err_l)];
    
    err_q = dudn_q_exact - dudn_q;
    err_q_L2 = [err_q_L2, sqrt(err_q' * err_q / length(err_q))];    % 离散l2误差
    % err_q_L2 = [err_q_L2,sqrt(err_q' * Mb_q * err_q)];            % 加权L2误差

    fprintf('\n=====h = %f=====\n',h(i));
    fprintf('err_l_L2 = %.4e,err_q_L2 = %.4e \n'...
        ,err_l_L2(i),err_q_L2(i));
    
    if i > 1
        order_l = log(err_l_L2(i-1) / err_l_L2(i)) / log(h(i-1) / h(i));
        order_q = log(err_q_L2(i-1) / err_q_L2(i)) / log(h(i-1) / h(i));
        fprintf('线性收敛阶 = %.2f,二次收敛阶 = %.2f \n',order_l,order_q);
    end
    % err_l_L2_pre = err_l_L2;
    % err_q_L2_pre = err_q_L2;
    
end
ord_l = polyfit(log(h),log(err_l_L2),1);
ord_q = polyfit(log(h),log(err_q_L2),1);
fprintf('\n 线性收敛阶 = %.2f\n',ord_l(1));
fprintf('二次收敛阶 = %.2f\n',ord_q(1));