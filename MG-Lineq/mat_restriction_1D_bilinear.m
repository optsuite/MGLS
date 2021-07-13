function R = mat_restriction_1D_bilinear(lev)
%  Restriction operator, coarse level dimension: 2^lev + 1

nx_l = 2^lev + 1;   nx_u = 2^(lev+1) + 1;

R = sparse(nx_l, nx_u);
R(1,1:2) = [2 1]/4;
R(end, end-1:end) = [1 2]/4;

mat = [1 2 1]/4;
dj = 2:4;
for di = 2:nx_l-1
    R(di, dj) = mat;
    dj = dj + 2;
end
