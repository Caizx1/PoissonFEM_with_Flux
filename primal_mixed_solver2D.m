function varargout = primal_mixed_solver2D(varargin)
%PRIMAL_MIXED_SOLVER2D  原始-混合格式求解 Poisson 方程
%   -Δu = f in Ω,  u = g on ∂Ω
%
%   新接口（结构体输出，nargout <= 1）:
%     result = primal_mixed_solver2D(geom, f, g, h_in, h_bd)
%     result = primal_mixed_solver2D(geom, f, g, h_in, h_bd, order)
%     result = primal_mixed_solver2D(geom, f, g, h_in, h_bd, order, refine_opts)
%     result = primal_mixed_solver2D(params_struct)
%     result 字段:
%       result.mesh      : p(2×np), t(n×nt), e(2×ne), np, nt
%       result.boundary  : e_b_nodes(2×nb), e_boundary(2×ne_bd), ne, xi_nodes(2×n_xi)
%       result.sol       : u(np×1), xi(n_xi×1)
%
%   旧接口（位置输出，nargout >= 7，向后兼容）:
%     [u, xi, p, t, e, e_b_nodes, e_boundary] = primal_mixed_solver2D(...)

    % 解析并校验输入参数
    params = init_solver_params(varargin{:});

    % 网格生成
    [p, t, e] = generate_meshes(params);

    % 矩阵组装
    [A, F, B, G, e_b_nodes, e_boundary] = assemble_system(p, t, e, params);

    % 求解线性系统
    [u, xi] = solve_saddle_point(A, F, B, G);

    % 计算 xi 自由度对应的节点坐标
    xi_nodes = compute_xi_nodes(e_b_nodes, e_boundary, params.order);

    % 构建输出结构体
    mesh_out     = struct('p', p, 't', t, 'e', e, ...
                          'np', size(p,2), 'nt', size(t,2));
    boundary_out = struct('e_b_nodes', e_b_nodes, 'e_boundary', e_boundary, ...
                          'ne', size(e_boundary,2), 'xi_nodes', xi_nodes);
    sol_out      = struct('u', u, 'xi', xi);

    % 根据 nargout 分发
    if nargout <= 1
        varargout = {struct('mesh', mesh_out, 'boundary', boundary_out, 'sol', sol_out)};
    else
        varargout = {u, xi, p, t, e, e_b_nodes, e_boundary};
    end
end

% ==================== 内部子函数 ====================

function [p, t, e] = generate_meshes(params)
% 生成内部网格，输出 p(2×np), t(n×nt), e(2×ne)

    [p0, ~, t0] = initmesh(params.geom, "Hmax", params.h_in);
    p0 = p0';            % np0 × 2
    t0 = t0(1:3, :)';    % nt0 × 3

    if params.refine_opts.use
        [p, e, t] = boundary_concentrated_mesh(p0, t0, params.h_in, ...
            params.geom, params.refine_opts.C, params.refine_opts.max_iter);
        % p: np×2, e: 2×ne, t: nt×3
    else
        p = p0;
        t = t0;
        T = auxstructure(t);
        e = T.bdEdge';    % 2 × ne
    end

    if strcmp(params.order, 'quadratic')
        [p, t] = linear2quadMesh(p, t);
        T = auxstructure(t(:, 1:3));
        e = T.bdEdge';
    end

    p = p';   % 2 × np
    t = t';   % n × nt  (n=3 linear, n=6 quadratic)
end

function [A, F, B, G, e_b_nodes, e_boundary] = assemble_system(p, t, e, params)
% 组装刚度矩阵 A、载荷向量 F、边界耦合矩阵 B 和右端项 G
% 输入 p(2×np), t(n×nt), e(2×ne)

    e_b_nodes  = boundaryDivide(params.geom, params.h_bd);
    ne_node    = size(e_b_nodes, 2);
    e_boundary = [1:ne_node; 2:ne_node, 1];

    if strcmp(params.order, 'linear')
        A = StiffnessAssembler2D(p, t, @(x,y) 1);
        F = LoadAssembler2D(p, t, params.f);
        [B, G] = BoundaryCoupling2D_independent(p, e, e_b_nodes, params.g);
    else
        A = StiffnessAssembler2D(p, t, @(x,y) 1, 'quadratic');
        F = LoadAssembler2D(p, t, params.f, 'quadratic');
        [B, G] = BoundaryCoupling2D_discontinuous_linear(p', t', ...
            e_b_nodes, e_boundary, params.g);
    end
end

function [u, xi] = solve_saddle_point(A, F, B, G)
% 求解鞍点系统 [A, B; B', 0] * [u; -xi] = [F; G]

    np     = size(A, 1);
    nb_dof = size(B, 2);
    K   = [A, B; B', sparse(nb_dof, nb_dof)];
    rhs = [F; G];
    sol = K \ rhs;

    u  = full(sol(1:np));
    xi = full(-sol(np + 1:end));
end

function xi_nodes = compute_xi_nodes(e_b_nodes, e_boundary, order)
% 计算 xi 自由度对应的物理坐标 (2 × n_xi)
%   linear:    边中点（常数边界元，每条边 1 个自由度）
%   quadratic: 高斯点（间断线性边界元，每条边 2 个自由度）

    ne = size(e_boundary, 2);

    if strcmp(order, 'linear')
        xi_nodes = zeros(2, ne);
        for i = 1:ne
            idxA = e_boundary(1, i);
            idxB = e_boundary(2, i);
            xi_nodes(:, i) = (e_b_nodes(:, idxA) + e_b_nodes(:, idxB)) / 2;
        end
    else
        t1 = 0.5 - 0.5 / sqrt(3);
        t2 = 0.5 + 0.5 / sqrt(3);
        xi_nodes = zeros(2, 2 * ne);
        for i = 1:ne
            idxA = e_boundary(1, i);
            idxB = e_boundary(2, i);
            pA = e_b_nodes(:, idxA);
            pB = e_b_nodes(:, idxB);
            xi_nodes(:, 2*(i-1)+1) = (1 - t1) * pA + t1 * pB;
            xi_nodes(:, 2*(i-1)+2) = (1 - t2) * pA + t2 * pB;
        end
    end
end
