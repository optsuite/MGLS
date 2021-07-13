function    u = Dirichlet_opt_control(u, pars)

%-------------------------------------------------
% Test problem 1
% boundary data
% u(x) = 0 on boundary

%-------------------------------------------------
% Test problem 2
% boundary data
% u0 = f(x)   , y = 0, 0 <= x <= 1
%      f(y)   , x = 0, 0 <= y <= 1
%      f(x)   , y = 1, 0 <= x <= 1
%      f(y)   , x = 1, 0 <= y <= 1
%   where f(x) = x(1-x) 


u(1,:) = 0;
u(end,:) = 0;
u(:,1) = 0;
u(:,end) = 0;