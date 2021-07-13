function H = elliptic_nonlinear_2D_3_hessian(u, pars)

%--------------------------------------------------------------------------
% nonlinear PDE
%
%   \Omega = [0,1]*[0,1]
%   - div grad u - \lambda u exp(u) = f
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

[nx, ny] = size(u);
dx = 1/(nx-1);  dy = 1/(ny-1);

area = 0.5*dx*dy;
% area = dx*dy;

lambda = 10;
%-------------------------------------------------
% set right hand side b
%   



lambda_expu = lambda * exp(u);
lambda_u_expu = lambda_expu.*u;

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


%     % Set dimension of hessian
%     nx2 = (nx - 2);   ny2 = (ny - 2);
% 
%     dx2 = 1/(dx)^2;     dy2 = 1/(dy)^2;     dxy = 1/dx/dy;
% 
%     % forward difference---------------------------------------------------
%     Hxx = dx2 * ones(nx-1,nx-1);
%     Hxy = dxy * zeros(nx-1,nx-1);
%     Hyy = dy2 * ones(nx-1,nx-1);
%     Hzz = lambda_u_expu(2:end-1,2:end-1) + lambda_expu(2:end-1,2:end-1);
% 
% 
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


