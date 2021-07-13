%FMG    Full Multigrid Solver
%
%       [X,RESIDS,ITS] = FMG(A,B) uses the full-multigrid cycle to 
%       solve the linear system AX=B.  RESIDS is the final residual
%       and ITS is 1.
%

function [x, rn,iter] = fmg(A,b, pars,rtol, max_it)

% build prolongation, restriction, and coefficient matricies
if ~isfield(pars,'have_pro_restr') || pars.have_pro_restr ~= 1
    pars = build_pro_restr_matrix( pars);
end

if ~isfield(pars,'have_coeff_mat') || pars.have_coeff_mat ~= 1
    pars = build_coeff_matrix(A, pars);
end


[x,rn,iter] = fmg_cycle(pars.lev_finest, b, pars, rtol, max_it);
