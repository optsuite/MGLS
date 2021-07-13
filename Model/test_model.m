function test_model

% Test the correctness of the gradient and the hessian by using automatic
% differentiation

% Required package: Intlab

% clc

nx = 2^3+1
% nx = 19
% rand('state', 0);


v = rand(nx); 
v = zeros(nx); 

% v(2,:) = v(2,:) + 1e-4;

% x = gradientinit(v);
x = hessianinit(v);
%---------------------------------------------------------------------
%---------------------------------------------------------------------
pars.model = 'minsurf_2D_1_fun2'; pars.model_hess = 'minsurf_2D_1_hessian';
% pars.model = 'elliptic_nonlinear_2D_1_fun'; pars.model_hess = 'elliptic_nonlinear_2D_1_hessian';
% pars.model = 'elliptic_nonlinear_2D_2_fun'; pars.model_hess = 'elliptic_nonlinear_2D_2_hessian';
% pars.model = 'bratu_2D_1_fun'; pars.model_hess = 'bratu_2D_1_hessian';
% pars.model = 'nonconvex_2D_2_fun'; pars.model_hess = 'nonconvex_2D_2_hessian';
% pars.model = 'possion2D_1';  pars.model_hess = 'possion2D_1_hessian';
% pars.model = 'Newton_fun';  pars.model_hess = 'Newton_hessian';
% pars.model = 'nonconvex_two_wells'; pars.model_hess = 'nonconvex_two_wells_hess';


tic; [f, g] = feval(pars.model, v, pars); toc
tic; y =feval(pars.model, x, pars); toc

fprintf('\nModel name: %s\n', pars.model);
fprintf('f - y.x: %e \n', f - y.x);

G = g(2:end-1,2:end-1);
g2 = reshape(y.dx,nx,nx);  G2 = g2(2:end-1,2:end-1);
 
% full(g)
% full(g2)

% full(G)
% full(G2)

fprintf('norm(g - g2): %e \n', norm(G-G2, 'fro'));

% % % Test Hessian
tic; H = feval(pars.model_hess, v, pars); toc
H2 = y.hx;
H2([1:nx nx+1:nx:nx^2-nx   2*nx:nx:nx^2-nx nx^2-nx+1:nx^2] , :) = [];
H2(:,[1:nx nx+1:nx:nx^2-nx   2*nx:nx:nx^2-nx nx^2-nx+1:nx^2] ) = [];

% size(H)
% size(H2)

% idx = 1:10; 
% full(H(idx,idx))
% full(H2(idx,idx))
fprintf('norm(H-H2): %e \n',norm(H - H2, 'fro'));

% eig(H)

% return;
%---------------------------------------------------------------------

% % Test the noncovex problem in Gratton, Sartenaer, Toint
% nx = 2^3+1
% % nx = 19
% % rand('state', 0);
% 
% gamma = rand(nx); u = rand(nx); v = [u, gamma];
% x = gradientinit(v);
% % x = hessianinit(v);
% pars.model = 'nonconvex_ls_fun_b';
% pars.model = 'nonconvex_ls_fun_9p';
% 
% tic; [f, g] = feval(pars.model, v, pars); toc
% tic; y =feval(pars.model, x, pars); toc
% 
% 
% fprintf('\nModel name: %s\n', pars.model);
% fprintf('f: %e, f - y.x: %e \n', f, f - y.x);
% 
%  G = [g{1}(2:end-1,2:end-1); g{2}(2:end-1,2:end-1)];
% 
%  g2 = reshape(y.dx,nx,2*nx);
%  gu = g2(:,1:nx); gr = g2(:,nx+1:end);
%  gu = gu(2:end-1,2:end-1);
%  gr = gr(2:end-1,2:end-1);  
%  G2 = [gu; gr];
%  
% % full(G)
% % full(G2)
% 
% fprintf('norm(g - g2): %e \n', norm(G-G2, 'fro'));
% % 
% 
% % Test Hessian
% % tic; [f, g, H] = feval(pars.model, v, pars); toc
% % tic; y =feval(pars.model, x, pars); toc
% % tic; H = feval(pars.model_hessian, v, pars); toc
% % H = y.hx;
% 
% % full(H)