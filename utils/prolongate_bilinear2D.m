
% Nine points prolongation
function u_p = prolongate_bilinear2D(u_l,pars)
    [nx ny] = size(u_l);    
    nx = 2*nx - 1;  ny = 2*ny-1; u_p = zeros(nx,ny);
    indx = [1:2:nx];  indy = [1:2:ny];  u_p(indx, indy) = u_l(1:end,1:end);
    indy = 2:2:ny-1;                    u_p(indx, indy) = 0.5*( u_p(indx, indy-1) + u_p(indx, indy+1) );
    indy = 1:ny; indx = 2:2:nx-1;       u_p(indx, indy) = 0.5*( u_p(indx-1, indy) + u_p(indx+1, indy) );
end