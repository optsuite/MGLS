function    u = Dirichlet_minsurf(u, pars)

%-------------------------------------------------
% Test problem 1
% boundary data
% u0 = f(x)   , y = 0, 0 <= x <= 1
%      0      , x = 0, 0 <= y <= 1
%      f(x)   , y = 1, 0 <= x <= 1
%      0      , x = 1, 0 <= y <= 1
%   where f(x) = x(1-x) 

%-------------------------------------------------
% Test problem 2
% boundary data
% u0 = f(x)   , y = 0, 0 <= x <= 1
%      f(y)   , x = 0, 0 <= y <= 1
%      f(x)   , y = 1, 0 <= x <= 1
%      f(y)   , x = 1, 0 <= y <= 1
%   where f(x) = x(1-x) 


% lev = pars.lev;
% nx = 2^lev + 1;     ny = nx;

nx = size(u,1);

dx = 1/(nx-1);
x = 0:dx:1;



if strcmp(pars.prob,'minsurf1') == 1
    vb = x.*(1-x);
 
    u(1,:) = vb;
    u(end,:) = vb;
    u(2:end-1, 1) = 0;
    u(2:end-1,end) = 0;
    
%     u = [vb; sparse(nx-2,ny); vb];

elseif strcmp(pars.prob,'minsurf2') == 1

    vb = x.*(1-x);

    u(1,:) = vb;
    u(end,:) = vb;
    u(2:end-1, 1) = vb(2:end-1)';
    u(2:end-1,end) = vb(2:end-1)';
        
%     u0 = [vb; vb(2:end-1)' sparse(nx-2,ny-2) vb(2:end-1)'; vb];


elseif strcmp(pars.prob,'minsurf3') == 1

    vb = 0.5*sin(2*pi*x);
%     vb = x.*(1-x);
%     u(1,:) = vb;
%     u(end,:) = vb;
%     u(2:end-1, 1) = -vb(2:end-1)';
%     u(2:end-1,end) = -vb(2:end-1)';

    u(1,:) = vb;
    u(end,:) = -vb;
    u(2:end-1, 1) = -vb(2:end-1)';
    u(2:end-1,end) = vb(2:end-1)';

%     u(1,:) = vb;
%     u(end,:) = vb;
%     u(2:end-1, 1) = vb(2:end-1)';
%     u(2:end-1,end) = vb(2:end-1)';    

%     u0 = [vb; -vb(2:end-1)' sparse(nx-2,ny-2) vb(2:end-1)'; -vb];
end
