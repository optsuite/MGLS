function    u = Dirichlet_nonconvex(u, pars)

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




if strcmp(pars.prob,'nonconvex_2D_1') == 1
    u(1,:) = 0;
    u(end,:) = 0;
    u(:,1) = 0;
    u(:,end) = 0;
elseif strcmp(pars.prob,'nonconvex_2D_2') == 1
    % u = 1 + x;
    nx = size(u,1);
    dx = 2/(nx-1);
    x = -1:dx:1;
    
    vb = 1+x;
    u(1,:) = 0;
    u(end,:) = 2;
    u(:,1) = vb;
    u(:,end) = vb;
    
elseif  strcmp(pars.prob,'nonconvex_2D_3') == 1
    % u = 1 - |x|
    nx = size(u,1);
    dx = 2/(nx-1);
    x = -1:dx:1;    
    vb = 1 - abs(x);
    %vb = 1 - x.^2;
%     vb = x.*(1-x.^2);

    u(1,:) = 0;
    u(end,:) = 0;
    u(:,1) = vb;
    u(:,end) = vb;
   
%     u(1,:) = vb;
%     u(end,:) = vb;
%     u(:,1) = 0;
%     u(:,end) = 0;    
end