function    u = Dirichlet_elliptic_1D_linear(u, pars)

%-------------------------------------------------
% Test problem 1
% boundary data



if pars.model_index == 1
%     vb = x.*(1-x);
 
    u(1,1) = 0;
    u(end,1) = 0;

   


elseif pars.model_index == 2

%     vb = x.*(1-x);

    u(1,1) = 0;
    u(end,1) = 0;
        

end
