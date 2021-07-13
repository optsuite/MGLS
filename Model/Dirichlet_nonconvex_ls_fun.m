function x = Dirichlet_nonconvex_ls_fun(x,pars)


% nx = size(x{1},1);
% dx = 1/(nx-1);
% x = 0:dx:1;
% vb = sin(2*pi*x);

x{1}(1,:) = 0;
x{1}(end,:) = 0;
x{1}(2:end-1, 1) = 0;
x{1}(2:end-1,end) = 0;

x{2}(1,:) = 0;
x{2}(end,:) = 0;
x{2}(2:end-1, 1) = 0;
x{2}(2:end-1,end) = 0;