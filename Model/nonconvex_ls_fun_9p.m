function [J, G] = nonconvex_ls_fun_9p(x, pars)

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
area = (dx*dy);

x = 0:dx:1;     y = 0:dy:1;
[X Y] = meshgrid(x,y);
u0 = sin(6*pi*X).*sin(2*pi*Y);
% u0 = cos(6*pi*X).*sin(2*pi*Y);

%-------------------------------------------------
uu0 = u - u0;

alpha = 1e3;

J = sum( sum( gamma.^2))/alpha ...
    + sum( sum( uu0.^2 ) );
 
vz = zeros(1,nx); vz1 = zeros(nx-1,1);

Jur = (4*([vz; u(1:nx-1,:)] + [u(2:nx,:); vz] + [vz', u(:,1:ny-1)] + [u(:,2:ny),vz']) + ...
        ([vz; vz1, u(1:nx-1,1:ny-1)] + [vz; u(1:nx-1,2:nx), vz1] + [vz1, u(2:nx,1:ny-1); vz] + [u(2:nx,2:nx), vz1; vz]) + ...
        -20*u )/(6*dx^2) - gamma.*u;    
J = J + sum( sum( Jur.^2 ) );
J = J*area;     
     
% compute gradient
if nargout > 1
    idx = 2:nx-1; idy = 2:ny-1;
    Gr = gamma(idx,idy)/(alpha/2) -2*(Jur(idx,idy).*u(idx,idy));
    Gu = 2*uu0(idx,idy) - 2*(Jur(idx,idy).*(gamma(idx,idy) + 10/(3*dx^2))) ...
         +4*(Jur(idx-1,idy) + Jur(idx+1,idy) + Jur(idx,idy-1)  + Jur(idx,idy+1) )/(3*dx^2) ...
         +(Jur(idx-1,idy-1) + Jur(idx+1,idy-1) + Jur(idx-1,idy+1)  + Jur(idx+1,idy+1) )/(3*dx^2) ;
    
    G{1,1} = area * [zeros(1,nx); zeros(nx-2,1) Gu zeros(nx-2,1); zeros(1,nx)];
    G{2,1} = area * [zeros(1,nx); zeros(nx-2,1) Gr zeros(nx-2,1); zeros(1,nx)];
end

