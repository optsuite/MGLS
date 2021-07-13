function R = mat_restriction_1D_bilinear_int(lev)
%  Restriction operator, coarse level dimension: 2^lev + 1

nx_l = 2^lev + 1;   nx_u = 2^(lev+1) + 1;

R = sparse(nx_l-2, nx_u-2);

mat = [1 2 1]/4;
dj = 1:3;
for di = 1:nx_l-2
    R(di, dj) = mat;
    dj = dj + 2;
end
