function [J, G, H] = elliptic_nonlinear_2D_3_fun(u, pars)

%--------------------------------------------------------------------------
% nonlinear PDE
%
%   \Omega = [0,1]*[0,1]
%   - div grad u + \lambda u exp(u) = f
%           u = 0 on boundary
%
%  ===> min J(u) = \int  0.5 *|grad u|^2  + \labmda (u exp(u) - exp(u) ) - f u dx
%   
%   gradient, hessian: only respect to interior points
%
%  source: 
%     Van Emden Henson;
%           Multigrid methods for nonlinear problems, an overview
%--------------------------------------------------------------------------

[nx, ny] = size(u); dx = 1/(nx-1);  dy = 1/(ny-1);

area = 0.5*dx*dy;
area = dx*dy;

lambda = 10;
%-------------------------------------------------
% set right hand side b
%   
x = 0:dx:1;     y = 0:dy:1; [X Y] = meshgrid(x,y);

% case 1: 
% exact solution:  (x - x^2)*(y - y^2);
% X2 = X - X.^2;      Y2 = Y - Y.^2;      XY2 = X2.*Y2; 
% f = 2*( X2 + Y2 ) + lambda*XY2.*exp(XY2);

% case 2: 
% exact solution:  (x^2 - x^3)*sin(3*pi*y);
X2 = (1-X).*X.^2;   Y2 = sin(3*pi*Y);
f = ( (9*pi^2 + lambda*exp(X2.*Y2)).*X2 + 6*X -2 ).*Y2;

delta_u_x_fwd = ( u(2:end,1:end-1) - u(1:end-1,1:end-1) )/dx;
delta_u_y_fwd = ( u(1:end-1,2:end) - u(1:end-1,1:end-1) )/dy;

% delta_u_x_bwd = ( u(2:end,2:end) - u(1:end-1,2:end) )/dx;
% delta_u_y_bwd = ( u(2:end,2:end) - u(2:end,1:end-1) )/dy;

lambda_expu = lambda * exp(u);
lambda_u_expu = lambda_expu.*u;

f_u = (lambda_u_expu - lambda_expu) - f.*u; 
% forward difference
J = 0.5 * area* sum( sum( delta_u_x_fwd.^2 + delta_u_y_fwd.^2)) ...
         + area*sum( sum( f_u(1:end-1,1:end-1) ) );
 
% % backward difference     
% J = J + 0.5 * area * sum( sum( (delta_u_x_bwd.^2 + delta_u_y_bwd.^2)...
%          + 2*f_u(2:end, 2:end) ) ) ;

% compute gradient
if nargout > 1

    % forward difference
    Gx = delta_u_x_fwd;
    Gy = delta_u_y_fwd;
    Gz =  lambda_u_expu(2:end-1,2:end-1) - f(2:end-1,2:end-1);

    G = (1/dx) *( Gx(1:end-1,2:end) -  Gx(2:end,2:end) ) ...
        + (1/dy) * ( Gy(2:end,1:end-1) -  Gy(2:end,2:end) ) ...
        + Gz ;

%     % backward difference
%     Gx = delta_u_x_bwd;
%     Gy = delta_u_y_bwd;
% %     Gz =  - 0.5*lambda_expu(2:end-1,2:end-1);
% 
%     G = G + (1/dx) *( Gx(1:end-1,1:end-1) -  Gx(2:end,1:end-1) ) ...
%         + (1/dy) * ( Gy(1:end-1,1:end-1) -  Gy(1:end-1,2:end) ) ...
%         + Gz;  
    
    G = area * [zeros(1,nx); zeros(nx-2,1) G zeros(nx-2,1); zeros(1,nx)];

end

if nargout > 2
    % Set dimension of hessian
%     nx2 = (nx - 2);   ny2 = (ny - 2);

    dx2 = 1/(dx)^2;     dy2 = 1/(dy)^2;     dxy = 1/dx/dy;

    % forward difference---------------------------------------------------
    Hxx = dx2 * ones(nx-1,nx-1);
    Hxy = dxy * zeros(nx-1,nx-1);
    Hyy = dy2 * ones(nx-1,nx-1);
    Hzz = lambda_u_expu(2:end-1,2:end-1) + lambda_expu(2:end-1,2:end-1);


    % diagonal element ((i,j), (i,j))
    H1 = ( Hxx(1:end-1,2:end) +  Hxx(2:end,2:end) ) + ...
         ( Hyy(2:end,1:end-1) +  Hyy(2:end,2:end) ) + ...
            2 *  Hxy(2:end,2:end) + Hzz;
    H1 = H1(:);

    %  element ((i,j), (i+1,j))
    H2 =  -Hxx(2:end,2:end) -  Hxy(2:end,2:end);
    H2(end,:) = 0;
    H2 = H2(:);

    %  element ((i+1,j), (i,j+1))
    H3 =  Hxy(2:end,2:end);
    H3(end,:) = 0;
    H3 = H3(:);

    %  element ((i,j), (i,j+1))
    H4 =  - Hxy(2:end,2:end) - Hyy(2:end,2:end);
    H4 = H4(:);
    
    H = area * spdiags([H4 [0;H3(1:end-1)] H2 H1 [0;H2(1:end-1)]  [zeros(nx-2,1);H3(1:end-nx+2)] [zeros(nx-2,1);H4(1:end-nx+2)]], ...
        [-(nx-2) -(nx-3) -1 0 1 (nx-3) (nx-2)] , (nx-2)^2 , (ny-2)^2);
    
%     % diagonal element ((i,j), (i,j))
%     H1 = ( Hxx(1:end-1,2:end) +  Hxx(2:end,2:end) ) + ...
%          ( Hyy(2:end,1:end-1) +  Hyy(2:end,2:end) ) + ...
%             2 *  Hxy(2:end,2:end) + Hzz + ...
%          ( Hxx(1:end-1,1:end-1) +  Hxx(2:end,1:end-1) ) + ...
%          ( Hyy(1:end-1,1:end-1) +  Hyy(1:end-1,2:end) ) + ...
%          2 *  Hxy(1:end-1,1:end-1) + Hzz;
%      
%     H1 = H1(:);
% 
%     %  element ((i,j), (i+1,j))
%     H2 =  -Hxx(2:end,2:end) -  Hxy(2:end,2:end) ...
%           -Hxx(2:end,1:end-1) -  Hxy(2:end,1:end-1);
%     H2(end,:) = 0;
%     H2 = H2(:);
% 
%     %  element ((i+1,j), (i,j+1))
%     H3 =  Hxy(2:end,2:end) +  Hxy(2:end,2:end);
%     H3(end,:) = 0;
%     H3 = H3(:);
% 
%     %  element ((i,j), (i,j+1))
%     H4 =  - Hxy(2:end,2:end) - Hyy(2:end,2:end) ...
%           - Hxy(1:end-1,2:end) - Hyy(1:end-1,2:end);
%     H4 = H4(:);
% 
%     H = area * spdiags([H4 [0;H3(1:end-1)] H2 H1 [0;H2(1:end-1)]  [zeros(nx-2,1);H3(1:end-nx+2)] [zeros(nx-2,1);H4(1:end-nx+2)]], ...
%         [-(nx-2) -(nx-3) -1 0 1 (nx-3) (nx-2)] , (nx-2)^2 , (ny-2)^2);
% 


%     % backward difference---------------------------------------------------
%     Hxx = dx2 * ones(nx-1,nx-1);
%     Hxy = dxy * zeros(nx-1,nx-1);
%     Hyy = dy2 * ones(nx-1,nx-1);
% 
%     % diagonal element ((i,j), (i,j))
%     H1 = ( Hxx(1:end-1,1:end-1) +  Hxx(2:end,1:end-1) ) + ...
%          ( Hyy(1:end-1,1:end-1) +  Hyy(1:end-1,2:end) ) + ...
%          2 *  Hxy(1:end-1,1:end-1) + Gz;
%     H1 = H1(:);
% 
%     %  element ((i,j), (i-1,j))
%     H2 =  -Hxx(2:end,1:end-1) -  Hxy(2:end,1:end-1);
%     H2(end,:) = 0;
%     H2 = H2(:);
% 
%     %  element ((i-1,j), (i,j-1))
%     H3 =  Hxy(2:end,2:end);
%     H3(end,:) = 0;
%     H3 = H3(:);
% 
%     %  element ((i,j), (i,j-1))
%     H4 =  - Hxy(1:end-1,2:end) - Hyy(1:end-1,2:end);
%     H4 = H4(:);
% 
% 
%     H = H + spdiags([H4 [0;H3(1:end-1)] H2 H1 [0;H2(1:end-1)]  [zeros(nx-2,1);H3(1:end-nx+2)] [zeros(nx-2,1);H4(1:end-nx+2)]], ...
%         [-(nx-2) -(nx-3) -1 0 1 (nx-3) (nx-2)] , (nx-2)^2 , (ny-2)^2);
% 
%     H = area * H;
%


end

