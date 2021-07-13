function [J, G] = nonconvex_ls_fun_b(x, pars)

%--------------------------------------------------------------------------
% nonlinear PDE
%  ===> min J(r,u) = \int [ r(x,y)^2 + (u(x,y)-u0(x,y))^2 + (laplace(u) -
%  r(x,y) u(x,y) )^2 dx dy
%   
%  u = u{1}
%
%   gradient, hessian: only respect to interior points
%
%--------------------------------------------------------------------------

%[nx, ny] = size(x); nx = nx/2; u = x(1:nx,:);   gamma = x(nx+1:end,:);
% [nx, ny] = size(x); ny = ny/2; u = x(:,1:ny);   gamma = x(:,ny+1:end);

u = x{1};   gamma = x{2}; [nx, ny] = size(u);

dx = 1/(nx-1);  dy = 1/(ny-1);
area = 0.5*(dx*dy);

x = 0:dx:1;     y = 0:dy:1;
[X Y] = meshgrid(x,y);
u0 = sin(6*pi*X).*sin(2*pi*Y);
% u0 = cos(6*pi*X).*sin(2*pi*Y);

%-------------------------------------------------
uu0 = u - u0;

alpha = 1e3;

J = sum( sum( gamma.^2))/alpha ...
    + sum( sum( uu0.^2 ) );
 
idx = 2:nx-1; idy = 2:ny-1;
Jur = (u(idx,idy-1) + u(idx,idy+1) + u(idx-1,idy) + u(idx+1,idy) -4*u(idx,idy) )/dx^2 ...
        - gamma(idx,idy).*u(idx,idy);
J = J + sum( sum( Jur.^2 ) );

idx = 1; idy = 2:ny-1;
Jur_u = (u(idx,idy-1) + u(idx+1,idy) + u(idx,idy+1) -4*u(idx,idy) )/dx^2 ...
        - gamma(idx,idy).*u(idx,idy);
    
idx = nx; idy = 2:ny-1;    
Jur_b = (u(idx,idy-1) + u(idx-1,idy) + u(idx,idy+1) -4*u(idx,idy) )/dx^2 ...
        - gamma(idx,idy).*u(idx,idy);
    
idx = 2:nx-1; idy = 1;    
Jur_l = (u(idx-1,idy) + u(idx,idy+1) + u(idx+1,idy) -4*u(idx,idy) )/dx^2 ...
        - gamma(idx,idy).*u(idx,idy);

idx = 2:nx-1; idy = ny;    
Jur_r = (u(idx-1,idy) + u(idx,idy-1) + u(idx+1,idy) -4*u(idx,idy) )/dx^2 ...
        - gamma(idx,idy).*u(idx,idy);

idx = 1; idy = 1;    
Jur11 = (u(idx,idy+1) + u(idx+1,idy) -4*u(idx,idy) )/dx^2 ...
        - gamma(idx,idy).*u(idx,idy);
    
idx = 1; idy = ny;    
Jur12 = (u(idx,idy-1) + u(idx+1,idy) -4*u(idx,idy) )/dx^2 ...
        - gamma(idx,idy).*u(idx,idy);

idx = nx; idy = 1;    
Jur21 = (u(idx-1,idy) + u(idx,idy+1) -4*u(idx,idy) )/dx^2 ...
        - gamma(idx,idy).*u(idx,idy);

idx = nx; idy = ny;    
Jur22 = (u(idx-1,idy) + u(idx,idy-1) -4*u(idx,idy) )/dx^2 ...
        - gamma(idx,idy).*u(idx,idy);
    
J = J + sum(Jur_u.^2) +  sum(Jur_b.^2) + sum(Jur_l.^2) + sum(Jur_r.^2) + ...
    Jur11^2 + Jur12^2 + Jur21^2 + Jur22^2;  
    
    
J = J*area;     
     
% compute gradient
if nargout > 1
    
    Gr = gamma(2:end-1,2:end-1)/(alpha/2) -2*(Jur.*u(2:end-1,2:end-1));
    Gu = 2*uu0(2:end-1,2:end-1) - 2*(Jur.*(gamma(2:end-1,2:end-1) + 4/(dx)^2));
    Gu(1:end-1,:) = Gu(1:end-1,:)+ 2*Jur(2:end,:)/(dx)^2;
    Gu(2:end,:) = Gu(2:end,:)+ 2*Jur(1:end-1,:)/(dx)^2;    
    Gu(:,1:end-1) = Gu(:,1:end-1)+ 2*Jur(:,2:end)/(dx)^2;
    Gu(:,2:end) = Gu(:,2:end)+ 2*Jur(:,1:end-1)/(dx)^2;    


    Gu(1,:) = Gu(1,:) + 2*Jur_u/(dx)^2;
    Gu(end,:) = Gu(end,:) + 2*Jur_b/(dx)^2;
    Gu(:,1) = Gu(:,1) + 2*Jur_l/(dx)^2;
    Gu(:,end) = Gu(:,end) + 2*Jur_r/(dx)^2;

    
%     % last column
%     Gr(:,end) = Gr(:,end)  -2*(Jur(:,end).*u(2:end-1,end-1));
%     Gu(:,end) = Gu(:,end) - 2*(Jur(:,end).*(gamma(2:end-1,end-1) + 4/(dx)^2));
%     Gu(:,end-1) = Gu(:,end-1)+ 2*Jur(:,end)/(dx)^2;
%     Gu(2:end,end) = Gu(2:end,end)+ 2*Jur(1:end-1,end)/(dx)^2;    
%     Gu(1:end-1,end) = Gu(1:end-1,end)+ 2*Jur(2:end,end)/(dx)^2;
% 
% %     %the last corner
%     Gr(end,end) = Gr(end,end) -2*(Jur(end,end).*u(end-1,end-1));    
%     Gu(end,end) = Gu(end,end) - 2*(Jur(end,end).*(gamma(end-1,end-1) + 4/(dx)^2));
%     Gu(end,end-1) = Gu(end,end-1)+ 2*Jur(end,end)/(dx)^2;
%     Gu(end-1,end) = Gu(end-1,end)+ 2*Jur(end,end)/(dx)^2;    

    
    G{1,1} = area * [zeros(1,nx); zeros(nx-2,1) Gu zeros(nx-2,1); zeros(1,nx)];
    G{2,1} = area * [zeros(1,nx); zeros(nx-2,1) Gr zeros(nx-2,1); zeros(1,nx)];
end

