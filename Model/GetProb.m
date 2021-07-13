function pars = GetProb(name,pars)

pars.FMLS_int_op = 'bilinear';
pars.FMLS_int_op = 'cubic';

pars.explict_Hessian  = 0;   % 1 - Have an explicit Hessian, 0 - no Hessian
pars.explict_Hessian  = 1;

pars.HessianVector_meth  = 0; % 1 - Have an explicit Hessian
pars.HessianVector_meth  = 2;

pars.domain = [0 1 0 1];
switch name
    case { 'minsurf1'}
        pars.eval_fun = 'minsurf_2D_1_fun2';
        pars.eval_hessian = 'minsurf_2D_1_hessian';
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_minsurf';
        %             load minsurf1_pde_toolbox;
        %             u_pde=tri2grid(p,t,u,x,x);
        %             clear b e g p t u

    case { 'minsurf2'}
        pars.eval_fun = 'minsurf_2D_1_fun2';
        pars.eval_hessian = 'minsurf_2D_1_hessian';
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_minsurf';
        %             load minsurf2_pde_toolbox;
        %             u_pde=tri2grid(p,t,u,x,x);
        %             clear b e g p t u

    case { 'minsurf3'}
        pars.eval_fun = 'minsurf_2D_1_fun2';
        pars.eval_hessian = 'minsurf_2D_1_hessian';
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_minsurf';
        %             load minsurf3_pde_toolbox;
        %             u_pde=tri2grid(p,t,u,x,x);
        %             clear b e g p t u

    case {'possion2D_1'}
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_possion';
        %u_pde = X.*Y.*(1-Y).*(1-X);
        pars.eval_fun = 'possion2D_1';
        pars.eval_hessian = 'possion2D_1_hessian';

    % This function has to be rewritten    
%     case {'possion2D_2'}
%         pars.HaveBoundary = 1;
%         pars.boundary_value = 'Dirichlet_possion';
%         u_pde = sin(2*pi*X).*sin(2*pi*Y);
%         pars.eval_fun = 'possion2D_2';
%         pars.eval_hessian = 'possion2D_2_hessian';

    case  {'ell_nonlinear1'}
        pars.eval_fun = 'elliptic_nonlinear_2D_1_fun';
        %     pars.eval_fun = 'elliptic_nonlinear_2D_1_fun2';
        pars.eval_hessian = 'elliptic_nonlinear_2D_1_hessian';
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_elliptic_nonlinear';
        %             load elliptic-nonlinear_pde_toolbox.mat;
        %             u_pde=tri2grid(p,t,u,x,x);
        %             clear b e g p t u

    case  {'ell_nonlinear2'}
        pars.eval_fun = 'elliptic_nonlinear_2D_2_fun';
        pars.eval_hessian = 'elliptic_nonlinear_2D_2_hessian';
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_elliptic_nonlinear';

        %             load elliptic-nonlinear_pde_toolbox2.mat;
        %             u_pde=tri2grid(p,t,u,x,x);
        %             clear b e g p t u

    case  {'ell_nonlinear3'}
        pars.eval_fun = 'elliptic_nonlinear_2D_3_fun';
        pars.eval_hessian = 'elliptic_nonlinear_2D_3_hessian';
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_elliptic_nonlinear';
        %             u_pde = (X.^2 - X.^3).*sin(3*pi*Y);

    case {'bratu1'}
        pars.eval_fun = 'bratu_2D_1_fun';
        pars.eval_hessian = 'bratu_2D_1_hessian';
        pars.eval_hessian_vec_prod = 'bratu_2D_1_hessian_vec_prod';
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_bratu';
        %             load bratu_pde_toolbox.mat;
        %             u_pde=tri2grid(p,t,u,x,x);
        %             clear b e g p t u

    case {'nonconvex_2D_1'}
        pars.eval_fun = 'nonconvex_2D_2_fun';
        pars.eval_hessian = 'nonconvex_2D_2_hessian';
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_nonconvex';
        pars.domain = [-1 1 -1 1];
    case {'nonconvex_2D_2'}
        pars.eval_fun = 'nonconvex_2D_2_fun';
        pars.eval_hessian = 'nonconvex_2D_2_hessian';
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_nonconvex';
        pars.domain = [-1 1 -1 1];

    case {'nonconvex_2D_3'}
        pars.eval_fun = 'nonconvex_2D_2_fun';
        pars.eval_hessian = 'nonconvex_2D_2_hessian';
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_nonconvex';
        pars.domain = [-1 1 -1 1];

        pars.kappa_x = 1e-5;
        pars.FMLS_int_op = 'bilinear';
        %             pars.FMLS_int_op = 'cubic';

    case {'nonconvex_two_wells'}
        pars.eval_fun = 'nonconvex_two_wells';
        pars.eval_hessian = 'nonconvex_two_wells_hess';
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_nonconvex_two_wells';
        pars.domain = [-1 1 -1 1];
        

        pars.kappa_x = 1e-5;
        pars.FMLS_int_op = 'bilinear';
        %             pars.FMLS_int_op = 'cubic';
        
    case {'nonconvex_ls'}
        pars.eval_fun = 'nonconvex_ls_fun';
        pars.eval_fun = 'nonconvex_ls_fun_b';
        pars.eval_fun = 'nonconvex_ls_fun_9p';        
        pars.HaveBoundary = [1;1];
        pars.boundary_value = 'Dirichlet_nonconvex_ls_fun';

        %             pars.kappa_x = 1e-5;
        %pars.FMLS_int_op = 'bilinear';
        %             pars.FMLS_int_op = 'cubic';

        pars.explict_Hessian  = 0;   % 1 - Have an explicit Hessian, 0 - no Hessian
        % pars.explict_Hessian  = 1;
        pars.HessianVector_meth  = 0; % 1 - Have an explicit Hessian
        % pars.HessianVector_meth  = 2;

        %             load elliptic-nonlinear_pde_toolbox.mat;
        %             u_pde=tri2grid(p,t,u,x,x);
        %             clear b e g p t u

%     case {'Newton'}
%         pars.eval_fun = 'Newton_fun';
%         pars.eval_hessian = 'Newton_hessian';
%         pars.HaveBoundary = 1;
%         pars.boundary_value = 'Dirichlet_Newton_fun';
%         M = 1; a = 4*M/3;
%         pars.domain = [-a a -a a];


        %             pars.kappa_x = 1e-5;
        %             pars.FMLS_int_op = 'bilinear';
        %             pars.FMLS_int_op = 'cubic';


    case {'opt_control1'}
        pars.eval_fun = 'opt_control_1_fun';
        %             pars.eval_fun = 'opt_control_1_fun_Bw';
        %pars.eval_hessian = 'elliptic_nonlinear_2D_1_hessian';
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_opt_control';
        %             load elliptic-nonlinear_pde_toolbox.mat;
        %             u_pde=tri2grid(p,t,u,x,x);
        %             clear b e g p t u

        pars.explict_Hessian  = 0;   % 1 - Have an explicit Hessian, 0 - no Hessian
        % pars.explict_Hessian  = 1;
        pars.HessianVector_meth  = 0; % 1 - Have an explicit Hessian
        % pars.HessianVector_meth  = 2;


    case {'opt_control2'}
        pars.eval_fun = 'opt_control_1_fun_Bw';
        %pars.eval_hessian = 'elliptic_nonlinear_2D_1_hessian';
        pars.HaveBoundary = 1;
        pars.boundary_value = 'Dirichlet_opt_control';
        %             load elliptic-nonlinear_pde_toolbox.mat;
        %             u_pde=tri2grid(p,t,u,x,x);
        %             clear b e g p t u

        pars.explict_Hessian  = 0;   % 1 - Have an explicit Hessian, 0 - no Hessian
        % pars.explict_Hessian  = 1;
        pars.HessianVector_meth  = 0; % 1 - Have an explicit Hessian
        % pars.HessianVector_meth  = 2;

end

end
