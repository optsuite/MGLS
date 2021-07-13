function pars = MGLS_Unc_pars_set(x0, pars)

% Options for mg_opt_unc() 
%
%--------------------------------------------------------------------------
% DESCRIPTION
%--------------------------------------------------------------------------
%
% pars = mg_opt_unc_pars_set(pars)
%
% If pars is empty upon input, pars will be returned containing the default options  
%
% Alternatively, if pars is passed with some fields already defined, those
% fields will be checked for errors, and the remaining fields will be added
% and initialized to their default values.
%
% Copyright (c) 2007.  Zaiwen Wen
%
%--------------------------------------------------------------------------


% Build basic information on the variables. These information is suitable
% for all levels


%% pars.IsCellX, pars.Xinfo and pars.HaveBoundary are three control variables, they are used to
%% define the operations, such as, x+a*y, norm, dot(x,y)
% is input an array or a cell? Reformulate x as a vector
pars.IsCellX = 0;
if iscell(x0);  
    pars.IsCellX = 1; 
    if ~isvector(x0); x0 = x0(:); end
    pars.nVarX = length(x0);
end

% par.Xinfo = 0, x is a vector
%           = 1, x is a matrix
if pars.IsCellX == 0
    pars.Xinfo = 0;
    if ~isvector(x0); 
        pars.Xinfo = 1; 
        pars = AddLowerBoundedNumberOption(pars, 'domain', [0 1 0 1], -inf, inf, 'initial solution');
    else
        pars = AddLowerBoundedNumberOption(pars, 'domain', [0 1], -inf, inf, 'initial solution');
    end
elseif pars.IsCellX == 1
    pars.Xinfo = zeros(pars.nVarX,1);
    for di = 1: pars.nVarX
        if ~isvector(x0{di}); 
            pars.Xinfo(di) = 1; 
            pars = AddLowerBoundedNumberOption(pars, 'domain', [0 1 0 1], -inf, inf, 'initial solution');
        else
            pars = AddLowerBoundedNumberOption(pars, 'domain', [0 1], -inf, inf, 'initial solution');
        end
    end
end

% Check boundary information
if ~isfield(pars, 'HaveBoundary')
   error('Please specify the type of boundary for each variable');
else
   if ~isfield(pars, 'boundary_value')
       error('Please give a function which provide the boundary condition');
   end
end

pars = AddLowerBoundedNumberOption(pars, 'explict_Hessian', 0, 0, 10, '1-explicit Hessian, 0-no Hessian, integer');
pars = AddLowerBoundedNumberOption(pars, 'HessianVector_meth', 0, 0, 10, '1-Hessian-vector product, 2-compute H*v by constructing H, 0-finite difference, integer');

if pars.explict_Hessian == 1 
    if (~isfield(pars, 'eval_hessian') || strcmp(pars.eval_hessian,'') == 1)
        error('no explicit Hessian, please set: pars.eval_hessian');
    end
elseif pars.explict_Hessian == 0
    if pars.HessianVector_meth == 1
        if (~isfield(pars, 'hessian_vec_prod') || strcmp(pars.hessian_vec_prod,'') == 1)
            fprintf('!!!!!no explicit Hessian, please set: pars.hessian_vec_prod !!!!!');
            pars.HessianVector_meth == 0;
        end
    elseif pars.HessianVector_meth == 2
        if (~isfield(pars, 'eval_hessian') || strcmp(pars.eval_hessian,'') == 1)
            %error('no explicit Hessian, please set: pars.eval_hessian');
            fprintf('!!!!!no explicit Hessian, please set: pars.hessian_vec_prod !!!!!');
            pars.HessianVector_meth == 0;
        end
    end
end


%-----------------------------------------------------------------------------
pars = AddLowerBoundedNumberOption(pars, 'lev_finest', 1, 1, 100, 'the coarsest level');
pars = AddLowerBoundedNumberOption(pars, 'lev_coarest', 1, 1, 100, 'the coarsest level');



pars = AddStringOption(pars, 'initial_x_coarse', 'old', {'old','restriction'}, 'initial point for coarse grid');


pars = AddLowerBoundedNumberOption(pars, 'display', 0, -10, 10, 'print information, -1=quiet, 0=some output, 1=more output. integer');
pars = AddLowerBoundedNumberOption(pars, 'fid', 1, -inf, inf, 'output file');


pars = AddStringOption(pars, 'enable_mg', 'yes', {'yes','no'}, 'do multigrid or not');
pars = AddLowerBoundedNumberOption(pars, 'maxiter', 1000, 1, 1e5, 'max number of iterations, integer');

if strcmp(pars.enable_mg, 'yes'); maxiter_coarse = 10; else   maxiter_coarse= 1000; end
pars = AddLowerBoundedNumberOption(pars, 'maxiter_coarse', maxiter_coarse, 1, 1e5, 'max number of iterations, integer');


pars = AddLowerBoundedNumberOption(pars, 'kappa_g', 1e-1, 0, 10, 'tolerance for switching MG step or direct step');
pars = AddLowerBoundedNumberOption(pars, 'kappa_g_a', 1e-5, 0, 10, 'tolerance for switching MG step or direct step');
pars = AddLowerBoundedNumberOption(pars, 'kappa_x', 1e-2, 0, 10, 'tolerance for two secutive x');
pars = AddLowerBoundedNumberOption(pars, 'eps_cvg', 1e-5, 0, 10, 'tolerance on the gradient');
pars = AddLowerBoundedNumberOption(pars, 'eps_cvg_min', 1e-16, 0, pars.eps_cvg, 'minimal tolerance on the gradient');
pars = AddLowerBoundedNumberOption(pars, 'f_rel_tol', 1e-14, 0, 10, 'tolerance on the relative change of the objective function');
pars = AddLowerBoundedNumberOption(pars, 'tolx', 1e-9, 0, 10, 'tolerance on the relative change of variables');

pars = AddLowerBoundedNumberOption(pars, 'eps_d', 1e-16, 0, 10, 'tolerance on the norm of the search direction');
pars = AddLowerBoundedNumberOption(pars, 'eps_stp', 1e-16, 0, 10, 'tolerance on the step size');
pars = AddLowerBoundedNumberOption(pars, 'good_stp', 0.8, 0.1, 1, 'good step size');


pars = AddLowerBoundedNumberOption(pars, 'iter_smooth', 1, 1, 100, 'number of smoothing');
pars = AddLowerBoundedNumberOption(pars, 'max_itr_direct', 5, 1, 100,'max number of itertations between two recursive steps, integer');
% set max_itr_direct to  8 if max iter on coarse is 10
% options for direct search


pars = AddStringOption(pars, 'direct_search', 'lbfgs', {'lbfgs','newton-fact','newton-pcg','newton-LMG'}, 'methods of direct search');
pars = AddStringOption(pars, 'ls_meth', 'BTMG', {'BTMG','Armijo-Wolfe'}, 'methods line search');

pars = AddLowerBoundedNumberOption(pars, 'pcg_mxitr', 20, 1, 1e5, 'max number of iterations of pcg, integer');
pars = AddLowerBoundedNumberOption(pars, 'pcg_tol', 1e-3, 0, 10, 'tolerance for pcg');
pars = AddLowerBoundedNumberOption(pars, 'tauh', 1e-6, 0, 10, 'perturbation for computing the H*x by finite difference');

pars = AddLowerBoundedNumberOption(pars, 'LMG_mxitr', 20, 1, 1e5, 'max number of iterations of linear multigrid method, integer');
pars = AddLowerBoundedNumberOption(pars, 'LMG_tol', 1e-3, 0, 1, 'tolerance for linear multigrid method');

pars = AddStringOption(pars, 'lbfgs_H_or_B', 'inverseB', {'inverseB','hessianB'}, 'L-BFGS method');
pars = AddLowerBoundedNumberOption(pars, 'lbfgs_inherit_coarse_info', 1, 0, 1, 'inherit information from the old minimization sequence');
pars = AddLowerBoundedNumberOption(pars, 'lbfgs_m', 5, 1, 100, 'number of storage for L-BFGS');
pars = AddLowerBoundedNumberOption(pars, 'eps_YS', 1e-20, 0, 10, 'tolerance on restarting L-BFGS, dot(Y,S)');
pars = AddLowerBoundedNumberOption(pars, 'eps_YS_r', 1e-8, 0, 10, 'tolerance on restarting L-BFGS, dot(Y,S)/(|Y|*|S|)');

% for hybrid method
pars = AddLowerBoundedNumberOption(pars, 'max_lbfgs_itr', 10, 1, 100, 'max number of iterations for L-BFGS in a sequel');

%----------------------------
% options for line search 
pars.pars_ls = [];
pars.pars_ls = AddLowerBoundedNumberOption(pars.pars_ls, 'maxiter', 100, 1, 1e5, 'max number of iterations of line search, integer');
pars.pars_ls = AddLowerBoundedNumberOption(pars.pars_ls, 'ftol', 1e-3, 0, 1, 'ftol for line search');
pars.pars_ls = AddLowerBoundedNumberOption(pars.pars_ls, 'gtol', 0.2, 0, 1, 'gtol for line search');
%----------------------------

pars.prolongation ='prolongation_regular';
pars.restriction ='restriction_regular';

pars = AddStringOption(pars, 'int_op_2D', 'prolongate_bilinear2D', {'prolongate_bilinear2D'}, '2D prolongation method');
pars = AddStringOption(pars, 'restrict_op_2D', 'restriction_bilinear2D', {'restriction_bilinear2D','restriction_simple2D'}, '2D restriction method');
pars = AddStringOption(pars, 'FMLS_int_op', 'bilinear', {'bilinear', 'cubic'}, '2D prolongation method for full approximation');

if pars.explict_Hessian  == 0 &&  strcmp(pars.direct_search, 'newton-fact') ==1
   error(' Newton method with factorization is required, but there is no explicit Hessian');
end

if pars.explict_Hessian  == 0 &&  strcmp(pars.direct_search, 'newton-LMG') ==1
   error(' Newton method with LMG is required, but there is no explicit Hessian');
end
