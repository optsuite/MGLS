function test_coarse_Hessian_quadratic_lineq

clc

lev = 5;
nx = 2^lev+1;
dx = 1/(nx-1);
% nx = 19
rand('state', 0);
xh = randn(nx);

%---------------------------------------------------------------------
pars.model = 'possion2D_1';

pars.lev_finest = lev;
pars.lev_coarest = 2;

% %-------------------------------------------------------------------------
% % full matrix
% 
tic; [fh, gh, Gh] = feval(pars.model,xh, pars); toc

b = rand(nx-2);    b = b(:);
x0 = zeros((nx-2)^2, 1);

% solve Gh * x = b
tol = 1e-3;  maxit = 50;
tic; [x,flag,relres,iter] = pcg(Gh,b,tol, maxit);toc
fprintf('PCG: flag: %d, relres: %3.2e, iter: %d\n\n', flag, relres, iter);

% profile on
tic; [x,relres,iter] = vmg(Gh,b, x0 ,pars, tol, maxit); toc
% profile report
fprintf('vmg: relres: %3.2e, iter: %d\n\n', relres, iter);

% profile on
tic; [x,relres,iter] = fmg(Gh,b, pars, tol, maxit); toc
% profile report
fprintf('fmg: relres: %3.2e, iter: %d\n', relres, iter);


