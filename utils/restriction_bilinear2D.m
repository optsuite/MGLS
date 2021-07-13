% Nine points prolongation
function u_l = restriction_bilinear2D(u_p,pars)
%     R = mat_restriction_2D_bilinear(pars.lev-1);
%     ull = R*u_p(:);
    
    [nx ny] = size(u_p);
    u_pt = zeros(nx+2, ny+2);     u_pt(2:end-1,2:end-1) = u_p;     clear u_p;
    nx = (nx+1)/2;  ny = (ny+1)/2;    u_l = zeros(nx,ny);
   
    indx = [1:nx];  indy = [1:ny];
    u_l(indx, indy) = 1/16 * u_pt(2*indx-1, 2*indy-1) + 1/8 * u_pt(2*indx-1, 2*indy) + 1/16 * u_pt(2*indx-1, 2*indy+1)...
                    + 1/8 * u_pt(2*indx, 2*indy-1)  + 1/4 * u_pt(2*indx, 2*indy) + 1/8 * u_pt(2*indx, 2*indy+1)...
                    + 1/16 * u_pt(2*indx+1, 2*indy-1) + 1/8 * u_pt(2*indx+1,2*indy) + 1/16 * u_pt(2*indx+1, 2*indy+1);
                
                
%     norm(ull-u_l(:))        
end
