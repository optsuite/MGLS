function [J, G] = nonconvex_ls_fun_1(x, pars)

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
[nx, ny] = size(x); ny = ny/2; u = x(:,1:ny);   gamma = x(:,ny+1:end);

%u = x{1};   gamma = x{2}; [nx, ny] = size(u);

dx = 1/(nx-1);  dy = 1/(ny-1);
area = 0.5*dx*dy;

x = 0:dx:1;     y = 0:dy:1;
[X Y] = meshgrid(x,y);
u0 = sin(6*pi*X).*sin(2*pi*Y);

%-------------------------------------------------
delta_u_x = ( u(2:end,:) - u(1:end-1,:) )/dx;
delta_u_y = ( u(:,2:end) - u(:,1:end-1) )/dy;

uu0 = u - u0;

J = sum( sum( gamma(2:end-1,2:end-1).^2))/1000 ...
    + sum( sum( uu0(2:end-1,2:end-1).^2 ) );
 
%J = sum( sum( gamma.^2))/1000 + sum( sum( (u - u0 ).^2 ) );
 
Lap_x = (delta_u_x(2:end,2:end-1) - delta_u_x(1:end-1,2:end-1))/dx;
Lap_y = (delta_u_y(2:end-1,2:end) - delta_u_y(2:end-1,1:end-1))/dy;
Jur = Lap_x + Lap_y - gamma(2:end-1,2:end-1).*u(2:end-1,2:end-1);
J = J + sum( sum( Jur.^2 ) ) ;
J = J*area;     
     
% compute gradient
if nargout > 1
    
    Gr = gamma(2:end-1,2:end-1)/500 -2*(Jur.*u(2:end-1,2:end-1));
    Gu = 2*uu0(2:end-1,2:end-1) - 2*(Jur.*(gamma(2:end-1,2:end-1) + 4/(dx)^2));
    Gu(1:end-1,:) = Gu(1:end-1,:)+ 2*Jur(2:end,:)/(dx)^2;
    Gu(2:end,:) = Gu(2:end,:)+ 2*Jur(1:end-1,:)/(dx)^2;    
    Gu(:,1:end-1) = Gu(:,1:end-1)+ 2*Jur(:,2:end)/(dx)^2;
    Gu(:,2:end) = Gu(:,2:end)+ 2*Jur(:,1:end-1)/(dx)^2;    

    
    G{1,1} = area * [zeros(1,nx); zeros(nx-2,1) Gu zeros(nx-2,1); zeros(1,nx)];
    G{2,1} = area * [zeros(1,nx); zeros(nx-2,1) Gr zeros(nx-2,1); zeros(1,nx)];
end

if nargout > 2
    % Set dimension of hessian
    nx2 = (nx - 2);   ny2 = (ny - 2);

    dx2 = 1/(dx)^2;     dy2 = 1/(dy)^2;     dxy = 1/dx/dy;

    % forward difference---------------------------------------------------
    Hxx = dx2 * ones(nx-1,nx-1);
    Hxy = dxy * zeros(nx-1,nx-1);
    Hyy = dy2 * ones(nx-1,nx-1);
    Hzz = lambda_u_expu(2:end-1,2:end-1) + lambda_expu(2:end-1,2:end-1);


    % diagonal element ((i,j), (i,j))
    H1 = ( Hxx(1:end-1,2:end) +  Hxx(2:end,2:end) ) + ...
         ( Hyy(2:end,1:end-1) +  Hyy(2:end,2:end) ) + ...
            2 *  Hxy(2:end,2:end) + Hzz + ...
         ( Hxx(1:end-1,1:end-1) +  Hxx(2:end,1:end-1) ) + ...
         ( Hyy(1:end-1,1:end-1) +  Hyy(1:end-1,2:end) ) + ...
         2 *  Hxy(1:end-1,1:end-1) + Hzz;
     
    H1 = H1(:);

    %  element ((i,j), (i+1,j))
    H2 =  -Hxx(2:end,2:end) -  Hxy(2:end,2:end) ...
          -Hxx(2:end,1:end-1) -  Hxy(2:end,1:end-1);
    H2(end,:) = 0;
    H2 = H2(:);

    %  element ((i+1,j), (i,j+1))
    H3 =  Hxy(2:end,2:end) +  Hxy(2:end,2:end);
    H3(end,:) = 0;
    H3 = H3(:);

    %  element ((i,j), (i,j+1))
    H4 =  - Hxy(2:end,2:end) - Hyy(2:end,2:end) ...
          - Hxy(1:end-1,2:end) - Hyy(1:end-1,2:end);
    H4 = H4(:);

    H = area * spdiags([H4 [0;H3(1:end-1)] H2 H1 [0;H2(1:end-1)]  [zeros(nx-2,1);H3(1:end-nx+2)] [zeros(nx-2,1);H4(1:end-nx+2)]], ...
        [-(nx-2) -(nx-3) -1 0 1 (nx-3) (nx-2)] , (nx-2)^2 , (ny-2)^2);




end

