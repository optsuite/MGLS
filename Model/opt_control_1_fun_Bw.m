function  [J, G] = opt_control_1_fun_Bw(u, pars)

%--------------------------------------------------------------------------
% optimal control
%
%   \Omega = [0,1]*[0,1]
%  min 0.5* || u - z ||^2 + 0.5*nu* || Bw||^2
%  s.t.   - div grad u - Bw = g
%           u = 0 on boundary
%  In this example, g = 0 and z = sin(2*pi*x) sin(pi*y)
%   B w = w , if in Omega_2 
%          = 0 , Omega\Omega_2
% 
%  ===> min J(u) = \int  0.5 *|| u - z ||^2  + 0.5*nu* || B(div grad u + g)||^2
%   
%   gradient only respect to interior points
%
%  source: 
%     A.Borzi, K.Kunisch;
%          A Multigrid Scheme for Elliptic Constrained Optimal Control
%          Problems
%--------------------------------------------------------------------------

[nx, ny] = size(u);
if nx~=ny
    error('Should be: nx == ny');
end

dx = 1/(nx-1);  dy = 1/(ny-1);

area = dx*dy;
nu = 1e-6;
%-------------------------------------------------
% set right hand side b
%   
x = 0:dx:1;     y = 0:dy:1;
[X Y] = meshgrid(x,y);

% case 1: 
% exact solution:  (x - x^2)*(y - y^2);
Z = sin(2*pi*X) .* sin(pi*Y);
LapU = zeros(nx,ny); 
LapU(2:end-1,2:end-1) = ( u(1:end-2, 2:end-1) + u( 3:end, 2:end-1)  ...
         +  u(2:end-1, 1:end-2) + u(2:end-1, 3:end)  ...
         -   4*u(2:end-1,2:end-1) )/dx^2;

W = (X - 0.5).^2 + (Y-0.5).^2 - sqrt(7/160)  < 0;   
LapU = W.*LapU;

% compute functional value     
J = 0.5*area* sum ( sum( ( u(2:end-1, 2:end-1) - Z(2:end-1, 2:end-1) ).^2 )) ;
J = J  + 0.5*nu*area* sum( sum( LapU.^2 ));

% compute gradient
if nargout > 1
    
    G = zeros(nx,ny);
    
    G(2:end-1, 2:end-1) =( u(2:end-1, 2:end-1) - Z(2:end-1, 2:end-1) ) - 4*nu*(LapU(2:end-1, 2:end-1) )/dx^2;
    
    G(2:end-1, 2:end-1) = G(2:end-1, 2:end-1) + ( LapU(1:end-2, 2:end-1) + LapU( 3:end, 2:end-1)  ...
         +  LapU(2:end-1, 1:end-2) + LapU(2:end-1, 3:end)  )* (nu/dx^2);
    
    G(2:end-1, 2:end-1) = G(2:end-1, 2:end-1)*area;
    
end


