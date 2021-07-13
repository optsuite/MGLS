function pars = build_coeff_matrix(A, pars)


% build restriction matrix:

pars.A_mat = cell(pars.lev_finest,1);
pars.L_mat = cell(pars.lev_finest,1);


pars.A_mat{pars.lev_finest,1} = A;
pars.L_mat{pars.lev_finest,1} = tril(A);
for di =pars.lev_finest-1:-1: pars.lev_coarest
    pars.A_mat{di,1} = pars.restriction_mat{di+1,1}*pars.A_mat{di+1,1}*pars.prolongation_mat{di+1, 1} ;
    pars.L_mat{di,1} = tril( pars.A_mat{di,1} ); ;
end






