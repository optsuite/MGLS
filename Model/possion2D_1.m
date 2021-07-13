function [J G H] = possion2D_1(u,pars)

%--------------------------------------------------------------------
% min J(u)
%   J(u) = 0.5* u'* A * u - u'*b 
%  
%   output: g and u has same dimensions
%
%   mesh: equispaced grid
%--------------------------------------------------------------------

% define parameters
[nx ny] = size(u); dx = 1/(nx-1); dy = 1/(ny-1);  area = dx*dy;

%-------------------------------------------------
% set right hand side b
%   f = 2(y(1-y) + x(1-x)) ==> u = xy(1-x)(1-y)
x = 0:dx:1;     y = 0:dy:1; [X Y] = meshgrid(x,y);
b = 2*(Y.*(1-Y) + X.*(1-X));    

delta_u_x_fwd = ( u(2:end,1:end-1) - u(1:end-1,1:end-1) )/dx;
delta_u_y_fwd = ( u(1:end-1,2:end) - u(1:end-1,1:end-1) )/dy;
% 
% delta_u_x_bwd = ( u(2:end,2:end) - u(1:end-1,2:end) )/dx;
% delta_u_y_bwd = ( u(2:end,2:end) - u(2:end,1:end-1) )/dy;



% forward difference
J =  0.5*sum( sum( delta_u_x_fwd.^2 + delta_u_y_fwd.^2 ) ) - sum( sum( u.*b ) );

J = area * J;

if nargout > 1
     
    Gx = delta_u_x_fwd;
    Gy = delta_u_y_fwd;
    Gz =  - b(2:end-1,2:end-1);

    G = (1/dx) *( Gx(1:end-1,2:end) -  Gx(2:end,2:end) ) ...
        + (1/dy) * ( Gy(2:end,1:end-1) -  Gy(2:end,2:end) ) ...
        + Gz ;

    G = area * [zeros(1,nx); zeros(nx-2,1) G zeros(nx-2,1); zeros(1,nx)];

end

if nargout > 2
     dx2 = 1/(dx)^2;     dy2 = 1/(dy)^2;    
     
    % forward difference---------------------------------------------------
    Hxx = dx2 * ones(nx-1,nx-1);
    Hyy = dy2 * ones(nx-1,nx-1);

    % diagonal element ((i,j), (i,j))
    H1 = ( Hxx(1:end-1,2:end) +  Hxx(2:end,2:end) ) + ...
         ( Hyy(2:end,1:end-1) +  Hyy(2:end,2:end) );
     
    H1 = H1(:);

    %  element ((i,j), (i+1,j))
    H2 =  -Hxx(2:end,2:end);
    H2(end,:) = 0;
    H2 = H2(:);

    %  element ((i,j), (i,j+1))
    H4 =  - Hyy(2:end,2:end);
    H4 = H4(:);

    H = area * spdiags([H4 H2 H1 [0;H2(1:end-1)]  [zeros(nx-2,1);H4(1:end-nx+2)]], ...
        [-(nx-2)  -1 0 1  (nx-2)] , (nx-2)^2 , (ny-2)^2);

    
    
end










