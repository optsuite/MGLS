function    u = Dirichlet_possion(u, pars)

%-------------------------------------------------
% Test problem 1
% boundary data


% lev = pars.lev;
% nx = 2^lev + 1;     ny = nx;

% nx = size(u,1);
% 
% dx = 1/(nx-1);
% x = 0:dx:1;



if strcmp(pars.prob,'possion2D_1') == 1
%     vb = x.*(1-x);
 
    u(1,:) = 0;
    u(end,:) = 0;
    u(2:end-1, 1) = 0;
    u(2:end-1,end) = 0;

elseif strcmp(pars.prob,'possion2D_2') == 1

%     vb = x.*(1-x);

    u(1,:) = 0;
    u(end,:) = 0;
    u(2:end-1, 1) = 0;
    u(2:end-1,end) = 0;
        

end
