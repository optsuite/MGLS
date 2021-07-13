function [J, G, H] = minsurf_2D_1_fun2(u, pars)

%--------------------------------------------------------------------------
% 2D minisurf
%
%  ===> min J(u) = \int sqrt( 1 +  (\nabla u)^2 )
%   
%   gradient, hessian: only respect to interior points
%--------------------------------------------------------------------------

[nx, ny] = size(u);
dx = 1/(nx-1);  dy = 1/(ny-1);

area = dx*dy;

%-------------------------------------------------
% set right hand side b
%   f = 2(y(1-y) + x(1-x)) ==> u = xy(1-x)(1-y)
% x = 0:dx:1;     y = 0:dy:1;
% [X Y] = meshgrid(x,y);


delta_u_x_fwd = ( u(2:end,1:end-1) - u(1:end-1,1:end-1) )/dx;
delta_u_y_fwd = ( u(1:end-1,2:end) - u(1:end-1,1:end-1) )/dy;
Valf_fwd = sqrt( 1 + delta_u_x_fwd.^2 + delta_u_y_fwd.^2 );
%Valf_fwd = sqrt( pars.minsurf_const + delta_u_x_fwd.^2 + delta_u_y_fwd.^2 );


% forward difference
J = area * sum( sum( Valf_fwd ) );
 
% % % backward difference  
% delta_u_x_bwd = ( u(2:end,2:end) - u(1:end-1,2:end) )/dx;
% delta_u_y_bwd = ( u(2:end,2:end) - u(2:end,1:end-1) )/dy;
% Valf_bwd = sqrt( 1 + delta_u_x_bwd.^2 + delta_u_y_bwd.^2 );

% J = J + area * sum( sum( Valf_bwd ));

% compute gradient
if nargout > 1

    % forward difference
    Valf_fwd = 1./Valf_fwd;
    Gx = delta_u_x_fwd .* Valf_fwd ;
    Gy = delta_u_y_fwd .* Valf_fwd;


    G = (1/dx) *( Gx(1:end-1,2:end) -  Gx(2:end,2:end) ) ...
        + (1/dy) * ( Gy(2:end,1:end-1) -  Gy(2:end,2:end) );

%     % backward difference
%     Valf_bwd = 1./Valf_bwd;
%     Gx = delta_u_x_bwd .* Valf_bwd ;
%     Gy = delta_u_y_bwd .* Valf_bwd;  
% 
%     G = G + (1/dx) *( Gx(1:end-1,1:end-1) -  Gx(2:end,1:end-1) ) ...
%         + (1/dy) * ( Gy(1:end-1,1:end-1) -  Gy(1:end-1,2:end) ) ;
    
    G = area * [zeros(1,nx); zeros(nx-2,1) G zeros(nx-2,1); zeros(1,nx)];
     
end

if nargout > 2
    dx2 = 1/(dx)^2;     dy2 = 1/(dy)^2;     dxy = 1/dx/dy;

    % forward difference---------------------------------------------------
    Valf_fwd2 = Valf_fwd.^3;
    Hxx = dx2* (Valf_fwd - Valf_fwd2.*(delta_u_x_fwd.^2));
    Hxy = dxy *( - Valf_fwd2.*delta_u_x_fwd.*delta_u_y_fwd);
    Hyy = dy2 * (Valf_fwd - Valf_fwd2.*(delta_u_y_fwd.^2) );
    
    
    % diagonal element ((i,j), (i,j))
    H1 = ( Hxx(1:end-1,2:end) +  Hxx(2:end,2:end) ) + ...
         ( Hyy(2:end,1:end-1) +  Hyy(2:end,2:end) ) + ...
         2 *  Hxy(2:end,2:end);
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

    H = spdiags([H4 [0;H3(1:end-1)] H2 H1 [0;H2(1:end-1)]  [zeros(nx-2,1);H3(1:end-nx+2)] [zeros(nx-2,1);H4(1:end-nx+2)]], ...
        [-(nx-2) -(nx-3) -1 0 1 (nx-3) (nx-2)] , (nx-2)^2 , (ny-2)^2);



%     % backward difference---------------------------------------------------
%     Valf_bwd2 = Valf_bwd.^3;
%     Hxx = dx2* (Valf_bwd - Valf_bwd2.*(delta_u_x_bwd.^2));
%     Hxy = dxy *( - Valf_bwd2.*delta_u_x_bwd.*delta_u_y_bwd);
%     Hyy = dy2 * (Valf_bwd - Valf_bwd2.*(delta_u_y_bwd.^2) );  
%     
%     % diagonal element ((i,j), (i,j))
%     H1 = ( Hxx(1:end-1,1:end-1) +  Hxx(2:end,1:end-1) ) + ...
%          ( Hyy(1:end-1,1:end-1) +  Hyy(1:end-1,2:end) ) + ...
%          2 *  Hxy(1:end-1,1:end-1);
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

    H = area * H;
%  
    

end

