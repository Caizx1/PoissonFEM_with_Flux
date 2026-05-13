function params = init_solver_params(varargin)
%INIT_SOLVER_PARAMS 解析并校验 primal_mixed_solver2D 的输入参数
%   params = init_solver_params(geom, f, g, h_in, h_bd)
%   params = init_solver_params(geom, f, g, h_in, h_bd, order)
%   params = init_solver_params(geom, f, g, h_in, h_bd, order, refine_opts)
%   params = init_solver_params(params_struct)
%
%   返回结构体，包含字段：
%       geom, f, g, h_in, h_bd, order, refine_opts(.use, .C, .max_iter)

    % ---- 调用约定分发 ----
    if nargin == 1 && isstruct(varargin{1})
        % 结构体语法
        opts = varargin{1};
        required = {'geom', 'f', 'g', 'h_in', 'h_bd'};
        for i = 1:length(required)
            if ~isfield(opts, required{i})
                error('init_solver_params:missingField', ...
                    '结构体缺少必要字段: ''%s''', required{i});
            end
        end
        params.geom = opts.geom;
        params.f    = opts.f;
        params.g    = opts.g;
        params.h_in = opts.h_in;
        params.h_bd = opts.h_bd;
        if isfield(opts, 'order') && ~isempty(opts.order)
            params.order = opts.order;
        end
        if isfield(opts, 'refine_opts')
            params.refine_opts = opts.refine_opts;
        end
    elseif nargin >= 5
        % 位置参数语法
        params.geom = varargin{1};
        params.f    = varargin{2};
        params.g    = varargin{3};
        params.h_in = varargin{4};
        params.h_bd = varargin{5};
        if nargin >= 6 && ~isempty(varargin{6})
            params.order = varargin{6};
        end
        if nargin >= 7
            params.refine_opts = varargin{7};
        end
    else
        error('init_solver_params:invalidArgs', ...
            '需要至少 5 个位置参数，或一个参数结构体。');
    end

    % ---- 默认值 ----
    if ~isfield(params, 'order') || isempty(params.order)
        params.order = 'linear';
    end
    if ~isfield(params, 'refine_opts')
        params.refine_opts = struct();
    end
    if ~isfield(params.refine_opts, 'use')
        params.refine_opts.use = false;
    end
    if ~isfield(params.refine_opts, 'C')
        params.refine_opts.C = 1.5;
    end
    if ~isfield(params.refine_opts, 'max_iter')
        params.refine_opts.max_iter = 10;
    end

    % ---- 校验 ----
    validateattributes(params.geom, {'numeric'}, {'nonempty'}, ...
        'init_solver_params', 'geom');
    validateattributes(params.f, {'function_handle'}, {}, ...
        'init_solver_params', 'f');
    validateattributes(params.g, {'function_handle'}, {}, ...
        'init_solver_params', 'g');
    validateattributes(params.h_in, {'numeric'}, {'positive', 'scalar'}, ...
        'init_solver_params', 'h_in');
    validateattributes(params.h_bd, {'numeric'}, {'positive', 'scalar'}, ...
        'init_solver_params', 'h_bd');
    params.order = validatestring(params.order, {'linear', 'quadratic'}, ...
        'init_solver_params', 'order');
end
