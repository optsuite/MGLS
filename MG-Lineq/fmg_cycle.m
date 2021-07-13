%FMG_CYCLE Full Multigrid Algorithm
%
%       U_OUT = FMG_CYCLE(LEVEL, B) uses the full-multigrid cycle to 
%       recursively solve the linear system AX=B at the given level.
%
%       No global variables are accessed.



function [u_out,rn,iter] = fmg_cycle(level, b, pars, rtol, max_it)

if level == pars.lev_coarest
   u_out   = pars.A_mat{level,1} \ b;
else 
   b_c	   = pars.restriction_mat{level,1}*b; 
   u_c_out = fmg_cycle(level-1, b_c, pars, rtol, max_it);
   u_f_in  = pars.prolongation_mat{level, 1}*u_c_out;
   
   pars.lev_finest = level;
        
   [u_out,rn,iter] = vmg(pars.A_mat{level,1},b,u_f_in,pars, rtol, max_it);
%    u_out   = vmg_cycle(level, b, u_f_in, pars);
end
