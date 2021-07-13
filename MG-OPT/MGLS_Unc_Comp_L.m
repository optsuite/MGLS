function [x, x0, output, stat] = MGLS_Unc_Comp_L(x0, g_up, pars, stat,  varargin)
%
% The purpose of "MGLS" is to demonstrate the algorithm proposed in 
%
%       @TECHREPORT{WenGoldfarb2007b,
%           author = {Wen, Zaiwen and Goldfarb, Donald},
%           title = {Line search Multigrid Methods for Large-Scale NonConvex Optimization},
%           institution = {Dept of IEOR, Columbia University},
%           year = {2007}
%       }
%
% Please refer to the example driver file "TestBatchMGDriver.m" and "GetProb.m" if the
% document is not clear.
%
%--------------------------------------------------------------------------
% GENERAL DESCRIPTION & INPUTS
%-------------------------------------------------------------------------- 
% 
% [x, x0, output, stat] = MGLS_Unc_Comp_L(x0, g_up, pars, stat,  varargin)
%
% is a recursive method to solve
%
% (P) min psi_l(x) = f_l(x) + (v_l)'*x, 
%   where {f_l(x_l)} is a hierarchy of discretizations of an infinite
%   dimensional minimization problem "min F(x)" on level "l". For example, the minimal
%   surface problem solves
%       min F(u) = \int_{Omega} sqrt( 1 +  (\nabla u)^2 ) dx
%       subject to, u = u_{Omega} on the boundary of Omega
%
% x0    initial solution. If x is a single variable in the continuous
%       problem F(x), then x is matrix or vector. If the continuous problem
%       has more than one variable, then x is a cell, where each field
%       contains one variable and the boundary condition of each variable
%       should be specified.
% 
% g_up  gradient from upper level. set g_up = [] at the beginning
%
% pars  a structure stores parameters to control the algorithm 
%  ------------------------------------------------------------------------
%  The following fields are required  
%
%   domain :  information of Omega
%   eval_fun : objective function value and gradient
%   HaveBoundary : scalar or vector, 0-without boundary, 1-Dirichlet boundary
%   boundary_value :  boundary value for each variable
%
%   lev_coarest : coarsest level 
%   lev_finest :  finest level 
%   enable_mg :   enable multigrid strategy or not, {'yes' or 'no'} 
%
%  ------------------------------------------------------------------------
%  The following fields are optional 
% 
%  'explict_Hessian': 1-explicit Hessian, 0-no Hessian, integer 
% 	default: 0, vaild range: [0, 10]
%
%   eval_hessian : the Hessian 
% 
% 	'HessianVector_meth': 1-Hessian-vector product, 2-compute H*v by constructing H, 0-finite difference, integer 
% 	default: 0, vaild range: [0, 10]
% 
% 	'lev_finest': the coarsest level 
% 	default: 1, vaild range: [1, 100]
% 
% 	'lev_coarest': the coarsest level 
% 	default: 1, vaild range: [1, 100]
% 
% 	'display': print information, -1=quiet, 0=some output, 1=more output. integer 
% 	default: 0, vaild range: [-10, 10]
% 
% 	'maxiter': max number of iterations, integer 
% 	default: 1000, vaild range: [1, 100000]
% 
% 	'maxiter_coarse': max number of iterations, integer 
% 	default: 10, vaild range: [1, 100000]
% 
% 	'kappa_x': tolerance for two secutive x 
% 	default: 0.1, vaild range: [0, 10]
% 
% 	'eps_cvg': tolerance on the gradient 
% 	default: 1e-05, vaild range: [0, 10]
% 
% 	'f_rel_tol': tolerance on the relative change of the objective function 
% 	default: 1e-16, vaild range: [0, 10]
% 
% 	'eps_d': tolerance on the norm of the search direction 
% 	default: 1e-16, vaild range: [0, 10]
% 
% 	'eps_stp': tolerance on the step size 
% 	default: 1e-16, vaild range: [0, 10]
% 
% 	'iter_smooth': number of smoothing 
% 	default: 1, vaild range: [1, 100]
% 
% 	'max_itr_direct': max number of itertations between two recursive steps, integer 
% 	default: 8, vaild range: [1, 100]
% 
% 	'pcg_mxitr': max number of iterations of pcg, integer 
% 	default: 20, vaild range: [1, 100000]
% 
% 	'pcg_tol': tolerance for pcg 
% 	default: 0.001, vaild range: [0, 10]
% 
% 	'tauh': perturbation for computing the H*x by finite difference 
% 	default: 1e-06, vaild range: [0, 10]
% 
% 	'LMG_mxitr': max number of iterations of linear multigrid method, integer 
% 	default: 20, vaild range: [1, 100000]
% 
% 	'LMG_tol': tolerance for linear multigrid method 
% 	default: 0.001, vaild range: [0, 10]
% 
% 	'lbfgs_inherit_coarse_info': inherit information from the old minimization sequence 
% 	default: 1, vaild range: [0, 1]
% 
% 	'lbfgs_m': number of storage for L-BFGS 
% 	default: 5, vaild range: [1, 100]
% 
% 	'max_itr_direct': max number of iterations for L-BFGS in a sequel 
% 	default: 10, vaild range: [1, 100]
% 
% 	'maxiter': max number of iterations of line search, integer 
% 	default: 100, vaild range: [1, 100000]
%
%   To set a field of opts, do  pars.[fieldname] = value. 
%
%
%--------------------------------------------------------------------------
% OUTPUTS
%--------------------------------------------------------------------------
% x             -exit solution
% x0            -store initial solution for next minimization sequence
% output        -a structure having the following fields
%    .exit      -exit status, 0='normal', 10='exceed max iterations'
%    .mesg      -a string describe the detailed exit status
%    .fval      -exit function value,
%    .norm_g    -exit l2 norm gradient
%    .pars      -options used
% stat          -a structure having the following fields
%    .nls       -num. of iterations on each level
%    .nfe       -num. of function evaluations on each level
%    .nge       -num. of gradient evaluations on each level
%    .nhe       -num. of hessian evaluations or matrix-vector products on each level
%    .n_cycle   -num. of multigrid cycles on each level
%
%--------------------------------------------------------------------------
% Copyright (c) 2008.  Zaiwen Wen 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Remark:
% *) Only a two-dimensional rectangular domain is supported in this version
% *) The condition of doing direct search can be either 
%       norm_g_low >= eps_cvg or norm_g_low >= eps_cvg_low
%    Sometimes, the latter is better
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Initialization
nx = sizeX(x0,pars);
if  pars.lev == pars.lev_finest
    % set a counter for the number of minimization sequences on each level
    if ~isfield(stat,'min_seq_counter') || isempty(stat.min_seq_counter)
        stat.min_seq_counter = zeros(pars.lev_finest,1);
        stat.min_seq_counter(pars.lev_finest,1) = 1;

        if strcmp(pars.initial_x_coarse, 'old') == 1
            stat.x_opt = cell(pars.lev_finest,1);
        end
    else
        stat.min_seq_counter = [stat.min_seq_counter;1];
        if strcmp(pars.initial_x_coarse, 'old') == 1
            stat.x_opt{pars.lev,1} = [];
        end
    end
else
    stat.min_seq_counter(pars.lev,1) = stat.min_seq_counter(pars.lev,1) + 1;
end


%--------------------------------------------------------------------------
% Initialization
%   Compute v
if pars.lev < pars.lev_finest && strcmp(pars.initial_x_coarse, 'old') == 1 ...
        && stat.min_seq_counter(pars.lev,1)-1 ~= 0
    %    stat.x_opt{pars.lev,1}
    if isempty(stat.x_opt{pars.lev,1})
        error('Initial x0 was not assigned');
    end
    x0 = stat.x_opt{pars.lev,1};
end
x = x0;

% compute function value and gradient
[f, g] = feval(pars.eval_fun, x, pars, varargin{:});
stat.nfe(pars.lev) = stat.nfe(pars.lev) + 1; stat.nge(pars.lev) = stat.nge(pars.lev) + 1;


if pars.lev == pars.lev_finest
    v_h = []; 
else
    % psi(x) = f(x) - v_h'*(x-x0)
    v_h = X_plus_Y(g, g_up, -1, pars); % v_h = g - g_up
    g = g_up;
end
f0 = f;     g0 = g;     norm_g = normG(g,'fro',pars);  

if pars.display >= 1
    fprintf(pars.fid, '%4s \t %4s \t %6s \t%1s \t %4s \t %6s \t %6s \t %10s \t %8s \t %7s %6s  \n', ...
        'Lev','Iter', 'd-Type', 'f', 'norm(g)', 'norm(d)', 'd''g','cos(theta)', 'stp', 'ls-Exit', 'ls_fun');
    fprintf(pars.fid, '%4d \t %4d \t %4s \t %+4.3e \t %4.3e \t %4.3e \t %+4.3e \t %4.3e \t %4.3e \t %5d \t %5d \n',...
        pars.lev, 0, ' ',  f, norm_g, [], [], [],[],[], []);
end

pars.step_method = 'MG'; stp = 1; d_old = []; sizeg = []; tauh = pars.tauh;
% Do not initialize istore after  LBFGS_init() since it will be assigned a value in LBFGS_init();
ndirect = 0; n_cycle = 0; x_cycle = [];   norm_x_cycle = inf;

% counter for updating LBFGS
if  strcmp(pars.direct_search, 'lbfgs') == 1
    istore = 0; status = 0; gamma = [];  ml = pars.lbfgs_m;
    LBFGS_init();
end

if pars.lev == pars.lev_finest; 
    vec_norm_g = zeros(pars.maxiter,1); x_type = zeros(pars.maxiter+2,1); vec_norm_g(1) = norm_g; 
    vec_f = zeros(pars.maxiter,1); vec_f(1) = f;
end

%--------------------------------------------------------------------------
% Begin main algorithm
% counter for the smoothing steps
iter_smooth = 1;

% set the max number of iterations for the current minimization sequence
% eps_cvg = pars.eps_cvg*5^( pars.lev_finest - pars.lev );
%eps_cvg = pars.eps_cvg*10^( pars.lev_finest - pars.lev );
eps_cvg = pars.eps_cvg*5^( -pars.lev_finest + pars.lev );
eps_cvg_low = pars.eps_cvg*5^( -pars.lev_finest + pars.lev -1 );

% eps_cvg = pars.eps_cvg*10^( -pars.lev_finest + pars.lev );
% eps_cvg = pars.eps_cvg;

% set tolerance, not smaller than eps_cvg_min
eps_cvg = max(eps_cvg, pars.eps_cvg_min);
iter = 0;
if norm_g < eps_cvg
    output.exit = 0; output.msg = 'first order optimal';
    % if converge as gradient norm, mark '**'
    pars.step_method = '**';         finalize_solu(); return
end

if pars.lev < pars.lev_finest
    maxiter = pars.maxiter_coarse;
else
    maxiter = pars.maxiter;
    %     maxiter = pars.maxiter *( pars.lev_finest - pars.lev + 1 );
end

for iter = 1:maxiter
    
    % assign f to fk
    fk = f; norm_gk = norm_g;
    stat.iter_count = stat.iter_count + 1;
    stat.iter_hist(stat.iter_count,1)  = pars.lev;

    if ( pars.lev > pars.lev_coarest ) && (iter_smooth <= pars.iter_smooth) && strcmp(pars.step_method,'DD') ~=1
        %mark smoothing step
        pars.step_method = 'SM';
        iter_smooth = iter_smooth + 1;
    else
        % reset the index of smoothing steps
        iter_smooth = 1;
        %******************************************************
        % Begin to compute search directionif any(pars.HaveBoundary); x0 = feval(pars.boundary_value, x0, pars); end
        if  (pars.lev > pars.lev_coarest) && strcmp(pars.enable_mg,'yes')==1 %&& iter > 1
            pars.var_type = 'dir'; g_low = feval(pars.restriction, g, pars);
            norm_g_low = normG(g_low, 'fro', pars);

            norm_x_x_last_cycle = inf;
            if n_cycle > 0
                norm_x_x_last_cycle = normX( X_plus_Y(x, x_cycle, -1, pars), 'fro', pars);
                if pars.display >= 2; fprintf('lev: %d, ||x-x_cycle||: %3.2e \n', pars.lev, norm_x_x_last_cycle); end
                %if pars.display >= 2; fprintf('lev: %d, ||x-x_cycle||: %3.2e, ||x_cycle||: %3.2e \n', pars.lev, norm_x_x_last_cycle, norm_x_cycle); end                
            end

            %if  ( ( norm_g_low >= pars.kappa_g * norm_g ) && (norm_g_low >= eps_cvg ))            
            %if  ( ( norm_g_low >= pars.kappa_g * norm_g ) && (norm_g_low >= eps_cvg )) && (norm_x_x_last_cycle >= pars.kappa_x || ndirect > pars.max_itr_direct ) ...    
            if  ( ( norm_g_low >= pars.kappa_g * norm_g ) && (norm_g_low >= eps_cvg )) && (norm_x_x_last_cycle >= pars.kappa_x *norm_x_cycle || ndirect > pars.max_itr_direct ) 
                pars.step_method = 'MG';
            else
                pars.step_method = 'DD';
            end

            if pars.display >= 2
                fprintf(pars.fid, 'd-Type: %4s \t norm_g: %4.3e \t norm_g_low: %4.3e \t  norm_g_low/norm_g: %4.3e \n', ...
                    pars.step_method,   norm_g, norm_g_low, norm_g_low/norm_g);
            end

        else
            pars.step_method = 'DD';
        end

    end

    % prepare for lbfgs
    if    strcmp(pars.direct_search, 'lbfgs') == 1
        % for Quasi-Newton's direction, try to accept the step size 1
        Y = g;
    end

    % compute a search direciton
    pars.norm_g = norm_g;
    if ( strcmp(pars.step_method, 'DD') == 1 ) || strcmp(pars.step_method,'SM')== 1
        ndirect = ndirect + 1;
        [gd, sizeg] = VecX(g, pars);
        switch pars.direct_search
            case {'newton-fact', 'newton-pcg', 'newton-LMG'}
                % for Newton's direction, try to accept the step size 1
                stp = 1.;
                [d, diter] = compute_direct_search_direction();
                % for iterative methods, we don't count the prepration of Hessian
                stat.nhe(pars.lev) = stat.nhe(pars.lev) + diter;
            case {'lbfgs'}
                % for Quasi-Newton's direction, try to accept the step size 1
                if status == 0; d = -gd;
                elseif status > 0; d = LBFGS_Hg_Loop(-gd); end
                %if status == 0; stp = 1e-3; else  stp = 1.;  end
                %if iter == 1; stp = 1/max(norm_g,1); else  stp = 1.;  end
                stp = 1.;
        end
        d = matX(d,sizeg,pars);

    else
        %  Compute the step recursively
        stp = 5.;        pars_low = pars;        pars_low.lev = pars.lev - 1;
        pars.var_type = 'var'; x_low = feval(pars.restriction, x, pars);
        if any(pars.HaveBoundary); x_low = feval(pars.boundary_value, x_low, pars); end

        % call multigrid recursively
        [x_low_opt, x_low, output_low, stat] = feval(pars.mg_opt,x_low, g_low,  pars_low, stat,  varargin{:});
        pars.var_type = 'dir'; d = feval(pars.prolongation, X_plus_Y(x_low_opt, x_low, -1, pars), pars);

        % record the point enter into the recursive step: direct-> direct-> recursive->new
        x_cycle = x;        n_cycle = n_cycle + 1;   x_type(iter+1) = 1;  ndirect = 0;
        norm_x_cycle = max(normX(x, 'fro', pars), 1);
        
        if pars.display >= 3
            gd_lo = ddotXY(g_low, X_plus_Y(x_low_opt,x_low, -1, pars), 1, pars);    gd_up = ddotXY(g, d, 1, pars);
            fprintf('gd_lo: %e \t gd_up: %e \t gd_up/gd_lo: %e \n', gd_lo, gd_up, gd_up/gd_lo);
        end
    end

    % if d is too short, it will not be able to make big progress, exit
    % issue: d maybe short at the very first step, it means the program
    %        return with the initial point. No recording in the figure of iteration
    %        iteration history.
    norm_d = normX(d, 'fro', pars);
    if ( norm_d < pars.eps_d )
        stat.lbfgs_seq_counter(pars.lev) = 0; % reset istore = 0;
        output.exit = 1; output.msg = 'd is too short';
        finalize_solu();  return;
    end

    % now, we have the search direction, begin to do line search
    switch pars.ls_meth
        case 'Armijo-Wolfe'
            % line search with the Armijo-Wolfe Condition, only for convex problem
            [stp, x, f, g, output_ls] = driver_ls_csrch_mg(stp, x, x0, f, g, d, v_h, pars, varargin{:});
        case 'BTMG'
            % line search using back tracking, inverse communication, for general non-convex problem
            % generate a stepsize which make the recursive direction is also a descent direction
            [stp, x, f, g, output_ls] = driver_ls_backtracking_mg(stp, x, x0, f, f0, g, g0, d, v_h, stat, pars, varargin{:});
            %[stp, x, f, g, output_ls] = driver_ls_BT_CSRCH_mg(stp, x, x0, f, f0, g, g0, d, v_h, stat, pars, varargin{:});            
    end

    % statistic for the line search
    output.ls_exit = output_ls.ls_exit;
    stat.nfe(pars.lev) = stat.nfe(pars.lev) + output_ls.nfe;
    stat.nge(pars.lev) = stat.nge(pars.lev) + output_ls.nge;

    % compute norm of gradient
    norm_g = normG(g, 'fro', pars);     %norm_g_inf = normG_v(g, 'inf', pars);
    if pars.lev == pars.lev_finest;  vec_f(iter+1) = f; vec_norm_g(iter+1) = norm_g; end

    d_deriv0 = output_ls.d_deriv0; theta = (abs(d_deriv0)/norm_gk/norm_d);
    % print information
    if pars.display >= 1
        fprintf(pars.fid, '%4d \t %4d \t %4s \t %+4.3e \t %4.3e \t %4.3e \t %+4.3e \t %4.3e \t %4.3e \t %5d \t %5d \n', ...
            pars.lev, iter, pars.step_method,  f, norm_g, norm_d, d_deriv0, theta, stp, output.ls_exit, output_ls.nfe);
    end

    % return if the step size is too small
    if ( stp < pars.eps_stp )
        stat.lbfgs_seq_counter(pars.lev) = 0; % reset istore = 0;
        output.exit = 1; output.msg = 'stp is too small';
        finalize_solu(); return;
    end
    
    % update storage for LBFGS
    if strcmp(pars.direct_search, 'lbfgs') == 1 %&& strcmp(pars.step_method,'MG') ~= 1
        if theta < 1e-2; % if d is a good direction
            istore = 0;
        end
        Save_LBFGS_Iter();
    end

    f_rel = (fk - f)/max( [abs(fk), abs(f), 1] );
    %fprintf('lev: %d, f_rel: %3.2e \n', pars.lev, f_rel);
    %******************************************************
    % Convergence test
    if norm_g < eps_cvg 
        output.exit = 0; output.msg = 'first order optimal';
        % if converge as gradient norm, mark '**'
        pars.step_method = '**';         finalize_solu(); return
    elseif  f_rel <  pars.f_rel_tol
        output.exit = 1; output.msg = 'no progress in objective';
        % if converge as gradient norm, mark '**'
        pars.step_method = '**';        finalize_solu(); return
    elseif stp*norm_d < pars.tolx && pars.lev == pars.lev_finest
        output.exit = 1; output.msg = 'no progress in x';
        % if converge as gradient norm, mark '**'
        pars.step_method = '**';        finalize_solu(); return
    end

end % end main algorithm

output.exit = 10; output.msg = 'Exceed Max MG Iterations'; finalize_solu();


% service routine
% prepare exit information
    function finalize_solu()
        output.fval = f;  output.norm_g = norm_g;
        stat.iter_count = stat.iter_count + 1;
        stat.iter_hist(stat.iter_count,1)  = pars.lev;
        stat.iter(pars.lev) = stat.iter(pars.lev) + iter;
        stat.n_cycle(pars.lev) = stat.n_cycle(pars.lev) + n_cycle;
        
        if pars.lev == pars.lev_finest;  
            stat.vec_norm_g = vec_norm_g(1:iter+1); stat.vec_f = vec_f(1:iter+1); stat.x_type = x_type(1:iter+1);
        end
        
        if strcmp(pars.initial_x_coarse, 'old') == 1
            stat.x_opt{pars.lev,1} = x;
        end
        
        if pars.display >= 2;  fprintf('lev: %d, ||x-x0||: %3.2e\n', pars.lev, normX( X_plus_Y(x,x0,-1,pars), 'fro', pars));  end
        if pars.display >= 1
            %fprintf(pars.fid, '%4d \t %4d \t %4s \t %+4.3e \t %4.3e** \t %4.3e \t %+4.3e \t %4.3e \t %4.3e \t %5d \t %5d \n', ...
            %    pars.lev, iter, pars.step_method,  f, norm_g, [], [], [], [], [],[]);
            fprintf(pars.fid, '%4d \t %4d \t %4s \t %+4.3e \t %4.3e* %s \n', ...
                pars.lev, iter, pars.step_method,  f, norm_g, output.msg);            
       end
    end

% second-order direction 
    function [d, diter] = compute_direct_search_direction()
        % Compute the step by Newton's method
        if pars.explict_Hessian  == 1
            H = feval(pars.eval_hessian, x, pars, varargin{:});    stat.nhe(pars.lev) = stat.nhe(pars.lev) + 1;
            switch pars.direct_search
                case {'newton-pcg'}
                    tol = pars.pcg_tol * norm_g;
                    [d,flag,relres,diter] = wpcg(H, -gd, tol, pars.pcg_mxitr, [],[]); d_old = d;
                    % R = cholinc(H,1e-3); [d,flag,relres,diter] = pcg(H, -gd, tol, pars.pcg_mxitr, R',R);
                    %     [d,flag,relres,iter] = pcg(@hessian_vec_prod_2D_second_order_pde, -g(:), tol, maxit);
                    %     H = feval(pars.eval_hessian_vec_prod, u, pars);
                    %     [d,flag,relres,iter] = pcg(@hessian_vec_prod_2D_second_order_pde, -g(:), tol, maxit, [], [], [], u, H);
                    %             [d,flag,relres,iter] = pcg(H, -g(:), tol, maxit);
                    diter = diter -1;
                case {'newton-LMG'}
                    if isempty(d_old); d_old = zeros(size(gd)); end
                    tol = pars.LMG_tol * norm_g;
                    lin_pars.lev_finest = pars.lev; lin_pars.lev_coarest = pars.lev_coarest;
                    lin_pars.display = pars.display;
                    lin_pars.restriction_mat = pars.restriction_mat;
                    lin_pars.prolongation_mat = pars.prolongation_mat;
                    pars.have_pro_restr = 1;
                    [d,relres,diter] = vmg(H, -gd, d_old,lin_pars, tol, pars.LMG_mxitr); d_old = d;
                    %         [d,relres,iter] = fmg(H, -g, lin_pars, tol, maxit);
                    if pars.display >= 1; fprintf('Lin-MG \t relres: %e \t iter: %d \n', relres, diter);  end
                    %fprintf('Lin-MG \t relres: %e \t iter: %d \n', relres, diter);
                case {'newton-fact'} %pars.lev == pars.lev_coarest ||
                    d = - H\gd;   diter = 0;
                    if gd'*d > 0
                        opts_eig.tol = 1e-4; opts_eig.disp = 0; opts_eig.issym = true; opts_eig.maxit = 100;
                        smin = eigs(H,1,'sa',opts_eig);
                        pp = 10; smin = 1.2*smin; if smin > 0; smin = -10; end
                        while pp > 0;  
                            fprintf('smallest eig value: %3.2e\n', smin);
                            [RR, pp] = chol(H-smin*speye(nx)); smin = 1.5*smin;  
                        end
                        d = -( RR\(RR'\gd));
                    end                    
%                     diter = 0;
%                     [RR, pp] = chol(H); 
%                     if pp == 0;
%                         d = -( RR\(RR'\gd));
%                     else
%                         opts_eig.tol = 1e-4; opts_eig.disp = 0; opts_eig.issym = true; opts_eig.maxit = 100;
%                         smin = eigs(H,1,'sa',opts_eig);
%                         %fprintf('smallest eig value: %3.2e\n', smin);
%                         smin = 1.2*smin; if smin > 0; smin = -10; end
%                         while pp > 0;  
%                             fprintf('smallest eig value: %3.2e\n', smin);
%                             [RR, pp] = chol(H-smin*speye(nx)); smin = 1.5*smin;  
%                         end
%                         d = -( RR\(RR'\gd));
%                     end
            end
        elseif pars.explict_Hessian  == 0
            %tol = pars.pcg_tol * norm_g;
            tol = 1e-8;
            if pars.HessianVector_meth  == 0;  
                cgHv = @HxFD;
            elseif pars.HessianVector_meth  == 1;
                cgHv = pars.hessian_vec_prod;
            elseif pars.HessianVector_meth  == 2;
                H = feval(pars.eval_hessian, x, pars, varargin{:});    cgHv = @(v) H*v;
            end
            [d, flag,relres, diter] = wpcg(cgHv, -gd, tol, pars.pcg_mxitr, [],[], []);
            if pars.display >= 1; fprintf('lev: %d, flag: %d, relres: %3.2e, diter: %d\n', pars.lev, flag, relres, diter); end
        end
    end

% approximate H*x by finite difference
    function y = HxFD(dh)
        if any(dh)
            dh = matX(dh,sizeg,pars); x_stp = X_plus_Y(x,dh, tauh,pars);
            if any(pars.HaveBoundary); x_stp = feval(pars.boundary_value, x_stp, pars); end
            [fh, gh] = feval(pars.eval_fun, x_stp, pars, varargin{:});
            stat.nfe(pars.lev) = stat.nfe(pars.lev) + 1; stat.nge(pars.lev) = stat.nge(pars.lev) + 1;
            if ~isempty(v_h);  gh = X_plus_Y(gh,v_h,-1,pars); end
            y = VecX(X_plus_Y(gh,g,-1,pars), pars)/tauh;
        else
            y = zeros(nx,1);
        end
    end

% initialize LBFGS
    function LBFGS_init()
        if  pars.lev == pars.lev_finest 
            istore = 0;
            if ~isfield(stat,'SK') || isempty(stat.SK)
                % no subfield SK, YK, STY, YTY
                stat.SK  = cell(pars.lev_finest,1);        stat.YK  = cell(pars.lev_finest,1);
                stat.perm = cell(pars.lev_finest,1);       stat.rho = cell(pars.lev_finest,1);  stat.gamma = cell(pars.lev_finest,1);
                % may get SK, YK, STY, YTY from last
                % set a counter for the number of storage in SK, YK, STY, YTY
                stat.lbfgs_seq_counter = zeros(pars.lev_finest,1);
            else
                % set a counter for the number of storage in SK, YK, STY, YTY
                stat.lbfgs_seq_counter = [stat.lbfgs_seq_counter; 0];
            end
            stat.SK{pars.lev,1}  = zeros(nx, pars.lbfgs_m);    stat.YK{pars.lev,1}  = zeros(nx, pars.lbfgs_m);
            stat.perm{pars.lev,1}  = [];      stat.rho{pars.lev,1}  = [];       stat.gamma{pars.lev,1}  = [];
        elseif pars.lev < pars.lev_finest 
            % If SK, YK, STY, YTY is empty at this level, create storage space for
            % them. Otherwise, use the old one. ( Currently, the objective function of different
            % minimization sequences on the coarse level are only different upto a
            % first order term, therefore, y_k = g_{k+1} - g_k should be coherent
            % for all objective functions on this level. )
            if pars.lbfgs_inherit_coarse_info == 1
                istore = stat.lbfgs_seq_counter(pars.lev);
            elseif pars.lbfgs_inherit_coarse_info == 0
                istore = 0;
            end
            if  istore == 0
                stat.SK{pars.lev,1}  = zeros(nx, pars.lbfgs_m);        stat.YK{pars.lev,1}  = zeros(nx, pars.lbfgs_m);
                stat.perm{pars.lev,1}  = [];        stat.rho{pars.lev,1}  = [];     stat.gamma{pars.lev,1}  = [];
            end
        end
        status = min(istore, pars.lbfgs_m);
    end

% computer y = H*v where H is L-BFGS matrix
    function y = LBFGS_Hg_Loop(dv)
        q = dv;   alpha = zeros(status,1);
        for di = status:-1:1;
            k = stat.perm{pars.lev,1}(di);
            alpha(di) = (q'*stat.SK{pars.lev,1}(:,k)) * stat.rho{pars.lev,1}(k);
            q = q - alpha(di)*stat.YK{pars.lev,1}(:,k);
        end
        y = stat.gamma{pars.lev,1}*q;
        for di = 1:status
            k = stat.perm{pars.lev,1}(di);
            beta = stat.rho{pars.lev,1}(k)* (y'* stat.YK{pars.lev,1}(:,k));
            y = y + stat.SK{pars.lev,1}(:,k)*(alpha(di)-beta);
        end
    end

% save storage for l-bfgs
    function Save_LBFGS_Iter()
        % after line search, g is updated to the gradient of the new point
        S = stp*VecX(d, pars); 
        Y = VecX(X_plus_Y(g, Y, -1, pars), pars); % Y = g - Y
        % save the computation of Y(:)'S(:)
        YS = Y'*S; norm_y = norm(Y); norm_s = norm(S); YS_ys = YS/(norm_y*norm_s);               
        if pars.display >= 2
            YS_yy = YS/(norm_y^2);  YS_ss = YS/(norm_s^2); 
            fprintf('lev: %d, YS_ys: %3.2e, YS_yy: %3.2e, YS_ss: %3.2e, YS: %3.2e, |Y|/|S|: %3.2e, status: %d\n', pars.lev, YS_ys, YS_yy,YS_ss, YS, norm_y/norm_s, status);
            %         Y = Y + ((2*(fk-f)+ GG'*S)/norm_s^2)*S;
            %         Y = Y + (1+max(-YS/norm_s^2, 0 ))*norm_g*S;
        end
        
        if YS < pars.eps_YS; istore = 0;  end % reset storage
        if  YS_ys > pars.eps_YS_r 
%         if YS > pars.eps_YS;
            istore = istore + 1;
            pos = mod(istore, ml); if pos == 0; pos = ml; end;
            stat.YK{pars.lev,1}(:,pos) = Y;          stat.SK{pars.lev,1}(:,pos) = S;
            stat.rho{pars.lev,1}(pos) = 1/(YS);      stat.gamma{pars.lev,1} = YS/norm_y^2;
            
            if istore <= ml; status = istore; stat.perm{pars.lev,1} = [stat.perm{pars.lev,1}, pos];
            else status = ml; stat.perm{pars.lev,1} = [stat.perm{pars.lev,1}(2:ml), stat.perm{pars.lev,1}(1)]; end
        
            stat.lbfgs_seq_counter(pars.lev) = min(istore, pars.lbfgs_m);
        end
    end
end