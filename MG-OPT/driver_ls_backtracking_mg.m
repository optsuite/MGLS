function [stp, x, f, g, output_ls] = driver_ls_backtracking_mg(stp, x, x0, f, f0, g, g0, d, v_h, stat, pars, varargin)

%--------------------------------------------------------------------------
% Backtracking line search algorithm
% inverse communication
%
% input:
%   stp     initial step size
%   x       crrent point
%   x0      the initial point of this minimizaiton sequence
%   f       current function value
%   f0      function value at x0
%   g       current gradient
%   g0      gradient at x0
%   d       search direction
%   v_h     first order addition to the objective function value
%   pars    parameters
%--------------------------------------------------------------------------

%reverse communication
workls.task = 1;
output_ls.iter = 0;    output_ls.nfe = 0;   output_ls.nge = 0;

%use Eucleadian norm
d_deriv0 = ddotXY(d, g, 1, pars); % d'*g
deriv_cd_g = 0;

output_ls.d_deriv0 = d_deriv0;
output_ls.f0 = f;
pars.pars_ls.maxiter = 10;
% pars.pars_ls.ftol = 1e-3;

% if pars.display >= 1;    pars.pars_ls.display = 'iter'; end
pars.pars_ls.display = 'no';

if pars.lev == pars.lev_finest
    while 1
        [stp, f, pars.pars_ls, workls] = ...
            ls_backtracking(stp, f, d_deriv0 , pars.pars_ls , workls);
        % Evaluate the function and the gradient at stp
        if (workls.task == 2)
            x_stp = X_plus_Y(x,d,stp,pars); % x_stp = x + stp*d
            if any(pars.HaveBoundary); x_stp = feval(pars.boundary_value, x_stp, pars); end

            f = feval(pars.eval_fun, x_stp, pars, varargin{:});
            % on finest level, f is the original function
            %f = f - feval(pars.inner_prod, v_h, x_stp - x0);
            output_ls.iter = output_ls.iter + 1;
            output_ls.nfe = output_ls.nfe + 1;
        else
            break;
        end
    end
    
   
elseif pars.lev < pars.lev_finest
    %     pars.pars_ls.ftol = 0.1;
    while 1
        [stp, f, pars.pars_ls, workls] = ...
            ls_backtracking_mg(stp, f, d_deriv0,  f0, deriv_cd_g, pars.pars_ls , workls);
        % Evaluate the function and the gradient at stp
        if (workls.task == 2)
            x_stp = X_plus_Y(x,d,stp,pars); % x_stp = x + stp*d
            if any(pars.HaveBoundary); x_stp = feval(pars.boundary_value, x_stp, pars); end
            f = feval(pars.eval_fun, x_stp, pars, varargin{:});
            xx0 = X_plus_Y(x_stp,x0,-1,pars); % x_stp - x0;
            f = f - ddotXY(v_h, xx0, 1, pars);
            %use Eucleadian norm
            deriv_cd_g = ddotXY(xx0, g0, 1, pars);
            %             ftol = pars.pars_ls.ftol;   gtol = 1 - ftol;
            %             fprintf('Level: %d, iter: %d, deriv: %5.4e \t dg0: %5.4e \t f: %5.4e \t stp: %5.4e \t fa: %5.4e \t fb: %5.4e\n', ...
            %                 pars.lev, output_ls.iter, deriv_cd_g, d_deriv0, f , stp, f1+stp*ftol*d_deriv0, f0 + gtol*deriv_cd_g);

            output_ls.iter = output_ls.iter + 1;
            output_ls.nfe = output_ls.nfe + 1;
        else
            break;
        end

    end

end

%     workls.msg
if workls.task ~= 0
    stp = 0;
    output_ls.d_deriv = d_deriv0;
    output_ls.ls_exit = workls.task;
else
    x = x_stp;
    % calculate gradient
    [f, g] = feval(pars.eval_fun, x_stp, pars, varargin{:});
    if ~isempty(v_h)
        f = f - ddotXY(v_h,  X_plus_Y(x_stp,x0,-1,pars), 1, pars);  % f = f - (v_h, x_stp - x0);
        g = X_plus_Y(g,v_h,-1,pars); %g - v_h;
    end
    output_ls.nge = output_ls.nge + 1;
end

% if pars.lev < pars.lev_finest && strcmp(pars.direct_search, 'lbfgs') == 1 && strcmp(pars.model_sec, 'lbfgsB')
%
% end

output_ls.d_deriv = d_deriv0;
output_ls.ls_exit = workls.task;


