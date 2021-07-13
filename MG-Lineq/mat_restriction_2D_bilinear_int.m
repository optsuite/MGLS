function R = mat_restriction_2D_bilinear_int(lev)
%  Restriction operator, coarse level dimension: 2^lev + 1

R = mat_restriction_1D_bilinear_int(lev);

R = kron(R, R);
