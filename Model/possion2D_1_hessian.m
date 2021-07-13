function H = possion2D_1_hessian(u,pars)



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





