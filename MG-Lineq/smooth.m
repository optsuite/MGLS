%SMOOTH Smooth a vector.
%
%       U_OUT = SMOOTH(LEVEL, B, U, FLAG) applies a smoother defined by the
%       global flag "smooth_flag" and the system AU=B to the vector U on the 
%       given grid level.  FLAG is set to 'pre', 'post', or 'coarse' and
%       defines the number of smoothings applied. 
%


function   u = smooth(level, b, u, flag, pars)


if strcmp(flag, 'pre') == 1
   nu = pars.nu1;
elseif  strcmp(flag, 'post') == 1
   nu = pars.nu2;
elseif strcmp(flag, 'coarse') == 1
   nu = 30;
end

   for i = 1:nu
       u = u + pars.L_mat{level,1}\(b - pars.A_mat{level,1}*u);
   end


