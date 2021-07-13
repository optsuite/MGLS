function [stp, x, f, g, output_ls] = driver_ls_csrch_mg(stp, x, x0, f, g, d, v_h, pars, varargin)

    %reverse communication
    workls.task = 1;
    output_ls.iter = 0;    output_ls.nfe = 0;
    d_deriv = ddotXY(d, g, 1, pars);
    output_ls.d_deriv0 = d_deriv;
    output_ls.f0 = f;
    
    while 1
        [stp, f, d_deriv, pars.pars_ls, workls] = ...
                ls_csrch(stp, f, d_deriv , pars.pars_ls , workls);
        % Evaluate the function and the gradient at stp
        if (workls.task == 2)
            x_stp = X_plus_Y(x,d,stp,pars); % x_stp = x + stp*d
            if any(pars.HaveBoundary); x_stp = feval(pars.boundary_value, x_stp, pars); end
            [f, g] = feval(pars.eval_fun, x_stp, pars, varargin{:});
            if ~isempty(v_h)
                f = f - ddotXY(v_h,  X_plus_Y(x_stp,x0,-1,pars), 1, pars);  % f = f - (v_h, x_stp - x0);
                g = X_plus_Y(g,v_h,-1,pars); %g - v_h;
            end
            d_deriv = ddotXY(d, g, 1, pars);
            output_ls.iter = output_ls.iter + 1;
            output_ls.nfe = output_ls.nfe + 1;
        else  % exit
            break
        end
    end
    x = x_stp;
    output_ls.d_deriv = d_deriv;
    output_ls.ls_exit = workls.task;            
    output_ls.nge = output_ls.nfe;
