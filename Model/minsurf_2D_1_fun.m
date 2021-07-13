function [J, G, H] = minsurf_2D_1_fun(u, pars)

%--------------------------------------------------------------------------
% 2D minisurf
%
%  ===> min J(u) = \int sqrt( 1 +  (\nabla u)^2 )
%   
%   gradient, hessian: only respect to interior points
%--------------------------------------------------------------------------

[nx, ny] = size(u);
dx = 1/(nx-1);  dy = 1/(ny-1);

% lev = log2(nx-1);
% hh = 10^(pars.lev_finest - lev+1);
hh = 1;

area = dx*dy;

delta_u_x = [( u(2,:) - u(1,:) )/(dx);
            ( u(3:end,:) - u(1:end-2,:) )/(2*dx);
            ( u(end,:) - u(end-1,:) )/(dx)];

delta_u_y = [ ( u(:,2) - u(:,1) )/(dy), ...
              ( u(:,3:end) - u(:,1:end-2) )/(2*dy), ...
              ( u(:,end) - u(:,end-1) )/(dy)];

Valf = sqrt( hh + delta_u_x.^2 + delta_u_y.^2 );


% forward difference
J = area*( sum( sum(Valf(2:end-1,2:end-1) ) ) + ...
    0.5*(sum( sum(Valf(1,2:end-1))) + sum(sum(Valf(end,2:end-1))) + ...
    sum( sum(Valf(2:end-1,1))) + sum(sum(Valf(2:end-1,end))) ) + ...
    0.25*(Valf(1,1) + Valf(1,end) + Valf(end,1) + Valf(end,end)) );

% compute gradient
if nargout > 1

    % forward difference
    Valf = 1./Valf;
    Gx = delta_u_x .* Valf ;
    Gy = delta_u_y .* Valf;

    G = zeros(nx-2,ny-2);
    % central FD
    G(2:end,:) = G(2:end,:) + (0.5/dx) *Gx(2:end-2,2:end-1);    
    G(1:end-1,:) = G(1:end-1,:) - (0.5/dx) *Gx(3:end-1,2:end-1);
    
    G(:,2:end) = G(:,2:end) + (0.5/dy) *Gy(2:end-1,2:end-2);        
    G(:,1:end-1) = G(:,1:end-1) - (0.5/dy) *Gy(2:end-1,3:end-1);    

    % (1,:)
    G(1,:) = G(1,:) + (0.5/dx) *Gx(1,2:end-1);
    G(end,:) = G(end,:) - (0.5/dx) *Gx(end,2:end-1);
    
    G(:,1) = G(:,1) + (0.5/dy) *Gy(2:end-1,1);
    G(:,end) = G(:,end) - (0.5/dy) *Gy(2:end-1,end);    
    
    G = area * [zeros(1,nx); zeros(nx-2,1) G zeros(nx-2,1); zeros(1,nx)];
     
end

if nargout > 2
    dx2 = 1/(dx)^2;     dy2 = 1/(dy)^2;     dxy = 1/dx/dy;

    % forward difference---------------------------------------------------
    Valf2 = Valf.^3;
    Hxx = dx2* (Valf - Valf2.*(delta_u_x.^2));
    Hxy = dxy *( - Valf2.*delta_u_x.*delta_u_y);
    Hyy = dy2 * (Valf - Valf2.*(delta_u_y.^2) );
    
    
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

