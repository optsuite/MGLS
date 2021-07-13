function    u = Dirichlet_nonconvex_two_wells(u, pars)


u(1,:) = 0;
u(end,:) = 0;
u(2:end-1, 1) = 0;
u(2:end-1,end) = 0;
    