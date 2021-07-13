%VMG_CYCLE V-Cycle algorithm.
%
%       U_OUT = VMG_CYCLE(LEVEL, B, U_IN) uses the V-cycle to recursively 
%       solve the linear system AX=B at the given level.  If the optional 
%       starting value U_IN is not passed then U_IN is set to 0's.
%


function u_out = vmg_cycle(level, b, u_in, pars)

% Use the zero vector for u_in as the default

% if nargin == 2,   
%    u_in = zeros(size(b));
% end


if level == pars.lev_coarest
   u_out   = pars.A_mat{level,1} \ b;
else 
   u       = smooth(level, b, u_in, 'pre',pars);
   r       = b - pars.A_mat{level,1}*u;
   
%    b_c     = pars.restriction_mat{level,1}*r;
%    u_c     = vmg_cycle(level-1, b_c , zeros(size(b_c)), pars);
%    correct = pars.prolongation_mat{level, 1}*u_c;
%    u       = u + correct;
%    norm_b_c = norm(b_c);
   
    norm_r = norm(r);

    coeff = 1e-2; 
    if norm_r > coeff * pars.rtol 
        b_c     = pars.restriction_mat{level,1}*r;
        norm_b_c = norm(b_c);
        if norm_b_c > (pars.rtol)*norm_r;
            u_c     = vmg_cycle(level-1, b_c , zeros(size(b_c)), pars);
            correct = pars.prolongation_mat{level, 1}*u_c;
            u       = u + correct;
        else
            [u_out,flag,relres,iter] = pcg(pars.A_mat{level,1}, b, pars.rtol, pars.max_it,[],[],u);
        end
        
    else % return
        %[u_out,flag,relres,iter] = pcg(pars.A_mat{level,1}, b, pars.rtol, maxit,[],[],u);
        
       u_out = u;
       return;
    end
   
   if (strcmp(pars.display, 'iter') == 1)
       fprintf('vmg-cycle \t level: %d \t norm_r: %e \t norm_b_c: %e \t tol: %e\n', ...
                                                level, norm_r, norm_b_c, pars.rtol);
   end
   
   u_out   = smooth(level, b, u, 'post',pars);
end
