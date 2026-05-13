close all;
clear;clc;

r = Rectg(0,0,1,1);
a = @(x,y)1;
f = @(x,y) 2 * pi^2 * sin(pi * x) .* sin(pi * y);
u_exact = @(x,y) sin(pi * x) .* sin(pi * y);

h = [0.2 0.1 0.05 0.025];
item = size(h,2);
err_l = zeros(item,1);
err_q = zeros(item,1);

for i = 1:item
    [pl,el,tl] = initmesh(r,"Hmax",h(i));
    [pq,tq] = linear2quadMesh(pl',tl(1:3,:)');
    pq = pq';tq = tq';

    Al = StiffnessAssembler2D(pl,tl,a);
    T = auxstructure(tl(1:3,:)');
    bdEdge = T.bdEdge;
    bdNodes = unique(bdEdge(:,1));
    nbdNodes = size(bdNodes,1);
    Al(:,bdNodes) = 0;Al(bdNodes,:) = 0;
    Al(bdNodes,bdNodes) = speye(nbdNodes);
    bl = LoadAssembler2D(pl,tl,f);
    bl(bdNodes) = 0;

    Aq = StiffnessAssembler2D(pq,tq,a,'quadratic');
    bq = LoadAssembler2D(pq,tq,f,'quadratic');

    [pqn,tqn] = trimeshRefine(pl',tl(1:3,:)');
    Tq = auxstructure(tqn);
    bdEdgeq = Tq.bdEdge;
    bdNodesq = unique(bdEdgeq(:));
    nbdNodesq = size(bdNodesq,1);

    Aq(bdNodesq,:) = 0; Aq(:,bdNodesq) = 0;
    Aq(bdNodesq,bdNodesq) = speye(nbdNodesq);
    bq(bdNodesq) = 0;

    ul = Al \ bl;
    err_l(i) = norm(u_exact(pl(1,:),pl(2,:))' - ul) / sqrt(size(ul,1));
    fprintf('使用线性元，h = %f，误差为%e \n',h(i),err_l(i));

    uq = Aq \ bq;
    err_q(i) = norm(u_exact(pq(1,:),pq(2,:))' - uq) / sqrt(size(uq,1));
    fprintf('使用二次元，h = %f，误差为%e \n',h(i),err_q(i));
end

ordl = polyfit(log(h),log(err_l),1);ordl = ordl(1);
ordq = polyfit(log(h),log(err_q),1);ordq = ordq(1);
fprintf('线性元收敛阶为%.4f，二次元收敛阶为%.4f \n',ordl,ordq);