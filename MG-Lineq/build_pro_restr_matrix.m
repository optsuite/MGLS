function pars = build_pro_restr_matrix(pars)


% build restriction matrix:
pars.restriction_mat  = cell(pars.lev_finest,1); 
pars.prolongation_mat = cell(pars.lev_finest,1); 

for di = pars.lev_coarest+1:pars.lev_finest
    pars.restriction_mat{di,1} =  mat_restriction_2D_bilinear_int(di-1);
    pars.prolongation_mat{di, 1} = 4*pars.restriction_mat{di,1}';
end






