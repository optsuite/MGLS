%VMG    Multigrid V-Cycle Solver
%
%       [X,RESIDS,ITS]=VMG(A,B,X0,RTOL,PRTOL,MAX_IT,MAX_TIME,MAX_MFLOP)
%       solves the system AX = B iteratively using multigrid cycles whose type
%       is defined by "cycle_flag".  The stopping criteria is given by
%       the input tolerances and limits.



%function [x,resids,its] = vmg(A,b,x,rtol,prtol,max_it,max_time,max_mflop)
function [x,rn,iter] = vmg(A,b,x,pars, rtol, max_it)

% build prolongation, restriction, and coefficient matricies
if ~isfield(pars,'have_pro_restr') || pars.have_pro_restr ~= 1
    pars = build_pro_restr_matrix( pars);
end

if ~isfield(pars,'have_coeff_mat') || pars.have_coeff_mat ~= 1
    pars = build_coeff_matrix(A, pars);
end

if ~isfield(pars,'display') 
    pars.display = 'no';
end


pars.nu1 = 1;
pars.nu2 = 1;

pars.rtol = rtol;
pars.max_it = max_it;

rn = norm(b - A * x);

%     c = (( rtol~=0 &  rn <= rtol*bn) ...
%        | (prtol~=0 & prn <= prtol*pbn) ...
%        | (iter >= max_it));

iter = 1;   
%while (~converged(bn,pbn,rn,prn,iter,rtol,prtol,max_it,max_time,max_mflop))
while  iter < max_it
    level = pars.lev_finest;
    x = vmg_cycle(level,b,x,pars);
    rn = norm(b - A*x);

    if (strcmp(pars.display, 'iter') == 1)
        fprintf('vmg \t norm_r: %e \t\n',  rn);
    end
    
    % convergence test
    if rn < rtol  
        return
    end
    
    iter = iter + 1;
end

