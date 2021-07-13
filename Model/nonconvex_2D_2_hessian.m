function  H =   nonconvex_2D_2_hessian(u, pars)

%--------------------------------------------------------------------------
% optimal control
%
%   \Omega = [-1,1]^2
%  min \int _{\Omega}  (1-ux^2)^2 + uy^2
%   s.t.  u = g  on boundary
%
%  If g = 0,  no minimizer
%  If g(x,y) = x + 1,   u = x + 1
%  If g(x,y) = 1 - |x|, u = 1 - |x|
% 
%   
%   gradient only respect to interior points
%
%  source: 
%       Rene Meziat
%       Analysis of two dimensional nonconvex variational problems
%--------------------------------------------------------------------------

[nx, ny] = size(u); 
if nx~=ny
    error('Should be: nx == ny');
end

dx = 2/(nx-1);  dy = 2/(ny-1);

area = dx*dy;
%-------------------------------------------------
% set right hand side b
%   
% x = -1:dx:1;     y = -1:dy:1;
% [X Y] = meshgrid(x,y);

delta_u_x_fwd = ( u(2:end,1:end-1) - u(1:end-1,1:end-1) )/dx;
%delta_u_y_fwd = ( u(1:end-1,2:end) - u(1:end-1,1:end-1) )/dy;

% delta_u_x_bwd = ( u(2:end,2:end) - u(1:end-1,2:end) )/dx;
% delta_u_y_bwd = ( u(2:end,2:end) - u(2:end,1:end-1) )/dy;



% Set dimension of hessian
dx2 = 1/(dx)^2;     dy2 = 1/(dy)^2;

% forward difference---------------------------------------------------
delta_u_x_fwd = delta_u_x_fwd.^2;
Hxx = (-4*dx2) *(1-delta_u_x_fwd) + (8*dx2) * delta_u_x_fwd;
Hyy = 2*dy2 * ones(nx-1,nx-1);

% diagonal element ((i,j), (i,j))
H1 = ( Hxx(1:end-1,2:end) +  Hxx(2:end,2:end) ) + ...
    ( Hyy(2:end,1:end-1) +  Hyy(2:end,2:end) );
H1 = H1(:);

%  element ((i,j), (i+1,j))
H2 =  -Hxx(2:end,2:end);
H2(end,:) = 0;
H2 = H2(:);

%  element ((i,j), (i,j+1))
H4 = -Hyy(2:end,2:end);
H4 = H4(:);

H = area * spdiags([H4 H2 H1 [0;H2(1:end-1)]  [zeros(nx-2,1);H4(1:end-nx+2)]], ...
    [-(nx-2) -1 0 1 (nx-2)] , (nx-2)^2 , (ny-2)^2);

% H = spdiags([H4 H2 H1 [0;H2(1:end-1)]  [zeros(nx-2,1);H4(1:end-nx+2)]], ...
%     [-(nx-2) -1 0 1 (nx-2)] , (nx-2)^2 , (ny-2)^2);


