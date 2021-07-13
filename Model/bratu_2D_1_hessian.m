function H = bratu_2D_1_hessian(u, pars)

%--------------------------------------------------------------------------
% 2D bratu equation:
%
%   \Omega = [0,1]*[0,1]
%   - div grad u - \lambda exp(u) = 0
%           u = 0 on boundary
%
%  ===> min J(u) = \int |grad u|^2  - \labmda exp(u) dx
%
%   gradient, hessian: only respect to interior points
%--------------------------------------------------------------------------

[nx, ny] = size(u);
dx = 1/(nx-1);  dy = 1/(ny-1);
% area = dx*dy;
area = 0.5*dx*dy;

%-------------------------------------------------
% set right hand side b
%   f = 2(y(1-y) + x(1-x)) ==> u = xy(1-x)(1-y)
% x = 0:dx:1;     y = 0:dy:1;
% [X Y] = meshgrid(x,y);

lambda = 6;


Gz =  - lambda * exp( u(2:end-1,2:end-1) );


nx2 = (nx - 2);   ny2 = (ny - 2);

dx2 = 1/(dx)^2;     dy2 = 1/(dy)^2;     dxy = 1/dx/dy;

% forward difference---------------------------------------------------
Hxx = dx2 * ones(nx-1,nx-1);
Hxy = dxy * zeros(nx-1,nx-1);
Hyy = dy2 * ones(nx-1,nx-1);


% diagonal element ((i,j), (i,j))
H1 = ( Hxx(1:end-1,2:end) +  Hxx(2:end,2:end) ) + ...
    ( Hyy(2:end,1:end-1) +  Hyy(2:end,2:end) ) + ...
    2 *  Hxy(2:end,2:end) + Gz + ...
    ( Hxx(1:end-1,1:end-1) +  Hxx(2:end,1:end-1) ) + ...
    ( Hyy(1:end-1,1:end-1) +  Hyy(1:end-1,2:end) ) + ...
    2 *  Hxy(1:end-1,1:end-1) + Gz;

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
H4 =  - Hxy(2:end,2:end) - Hyy(2:end,2:end) - Hxy(1:end-1,2:end) - Hyy(1:end-1,2:end);
H4 = H4(:);

H = area * spdiags([H4 [0;H3(1:end-1)] H2 H1 [0;H2(1:end-1)]  [zeros(nx-2,1);H3(1:end-nx+2)] [zeros(nx-2,1);H4(1:end-nx+2)]], ...
    [-(nx-2) -(nx-3) -1 0 1 (nx-3) (nx-2)] , (nx-2)^2 , (ny-2)^2);


