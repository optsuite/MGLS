function [x_mg, output, stat] = Full_MGLS_Unc(x0, pars,  varargin)
% The dirver file for the full multigrid method
%
%--------------------------------------------------------------------------
% Copyright (c) 2008.  Zaiwen Wen 
%--------------------------------------------------------------------------


if isempty(x0)
    error('please provide an initial solution!');
end

pars = MGLS_Unc_pars_set(x0, pars);

% if possible, build projection and prolongation matrix
if strcmp(pars.direct_search, 'newton-LMG') == 1
    pars = build_pro_restr_matrix( pars);
    pars.have_pro_restr = 1;
end

eps_cvg = pars.eps_cvg;
domain = pars.domain;
lev_coarest = pars.lev_coarest;
lev_finest = pars.lev_finest;

pars.lev = pars.lev_finest;

maxiter = pars.maxiter;
maxiter_coarse = pars.maxiter_coarse;
% pars.maxiter = maxiter_coarse;

stat.nfe = zeros(pars.lev,1);
stat.nge = zeros(pars.lev,1);
stat.nhe = zeros(pars.lev,1);
stat.norm_g = zeros(pars.lev,1);
stat.cpu = zeros(pars.lev,1);
stat.n_cycle = zeros(pars.lev,1);
stat.iter = zeros(pars.lev,1);
stat.iter_count = 1;
stat.iter_hist = [pars.lev_coarest];

stat.iter_gf = zeros(pars.lev,1);
stat.vec_norm_g = [];
stat.vec_f = [];
stat.x_type = [];

for iter = lev_coarest : lev_finest

    if pars.display >= 1
        fprintf(pars.fid, '\n \nFull Multigrid Method, Enter %d th \n\n', iter);
    end

%     pars.eps_cvg = eps_cvg*5^( -lev_finest + iter );
    
    pars.lev = iter;
    pars.lev_finest = pars.lev;
    pars.lev_coarest = lev_coarest;

    
    %     lev = iter + 1;   nx = 2^lev + 1;     ny = nx;
    [nx, ny] = size(x0);

%     if iter == lev_coarest
%         pars.maxiter = max(maxiter_coarse*5, floor(0.1*maxiter));
%         pars.maxiter = maxiter;
%         pars.maxiter_coarse = maxiter_coarse;
%     elseif iter == lev_finest
%         pars.maxiter = maxiter;
%         pars.maxiter_coarse = maxiter_coarse;
%     else
%         pars.maxiter = max(maxiter_coarse*5, floor(0.1*maxiter));
%         %         pars.maxiter = maxiter;
%         %         pars.maxiter = maxiter_coarse*5;
%     end

    %     %----------------------------------------------------------------------
    %     % plot initial solution
    %     dx = 1/(nx-1);
    %     x = 0:dx:1;
    %
    %     [X Y] = meshgrid(x,x);
    %     perf_fig = figure(iter+8);
    %     mesh(X, Y, x0)
    %     title(strcat('Initial Solution', ' level = ' , num2str(pars.lev)) )
    %     print(perf_fig , '-depsc2', strcat('.\results\',pars.eval_fun,'_No_',...
    %         num2str(pars.model_index),'_D_',pars.direct_highlevel,...
    %         '_SM_',num2str(pars.iter_smooth),'_mg_lev',num2str(pars.lev_finest),'.eps'));
    %     %----------------------------------------------------------------------

    % solve level p+1 by multigrid method


    stat1.nfe = zeros(pars.lev,1);
    stat1.nge = zeros(pars.lev,1);
    stat1.nhe = zeros(pars.lev,1);
    stat1.n_cycle = zeros(pars.lev,1);
    stat1.iter = zeros(pars.lev,1);
    stat1.iter_count = 0;
    stat1.iter_hist = zeros(10^3,1);

    %     pars.eps_cvg = pars.eps_cvg*5^( -pars.lev_finest + pars.lev );

    solve_t = cputime;
    [x_mg, x0, output, stat1] = feval(pars.mg_opt, x0, [], pars, stat1, varargin{:});
    stat.cpu(pars.lev) = cputime-solve_t;
    stat.norm_g(pars.lev) = output.norm_g;
    
    stat.nfe(1:pars.lev) = stat.nfe(1:pars.lev) + stat1.nfe(1:pars.lev);
    stat.nge(1:pars.lev) = stat.nge(1:pars.lev) + stat1.nge(1:pars.lev);
    stat.nhe(1:pars.lev) = stat.nhe(1:pars.lev) + stat1.nhe(1:pars.lev);
    stat.n_cycle(1:pars.lev) = stat.n_cycle(1:pars.lev) + stat1.n_cycle(1:pars.lev);
    stat.iter(1:pars.lev) = stat.iter(1:pars.lev) + stat1.iter(1:pars.lev);
    stat.iter_count = stat.iter_count + stat1.iter_count;
    stat.iter_hist = [stat.iter_hist; stat1.iter_hist(1:stat1.iter_count,1)];
    %     stat.vec_norm_g = stat1.vec_norm_g;
    %     stat.vec_f = stat1.vec_f;
    %     stat.x_type = stat1.x_type;

    if iter > lev_coarest
        stat.cpu(pars.lev) = stat.cpu(pars.lev) + stat.cpu(pars.lev-1);
        stat.iter_gf(pars.lev) = stat.iter_gf(pars.lev-1) + length(stat1.vec_norm_g);
        stat.vec_norm_g = [stat.vec_norm_g; stat1.vec_norm_g];
        stat.vec_f = [stat.vec_f; stat1.vec_f];
        stat.x_type = [stat.x_type; stat1.x_type];
    end
    
    %----------------------------------------------------------------------
    % plot difference of the solution and the initial points
    %     [nx, ny] = size(x_mg);
    %     dx = 1/(nx-1);
    %     x = 0:dx:1;
    %
    %     [X Y] = meshgrid(x,x);
    %     perf_fig = figure(iter+30);
    %     mesh(X, Y, x_mg - x0)
    %     title('x - x0');
    %----------------------------------------------------------------------

    if iter < lev_finest
        x0 = x_mg;

        % prolong initial point from level p to level p+1
        switch pars.FMLS_int_op
            case 'bilinear'
                x0 = feval(pars.prolongation, x0, pars);
            case 'cubic'
                nx = 2^iter + 1;
                dx = (domain(2)-domain(1))/(nx-1); x = domain(1):dx:domain(2);
                [X1 Y1] = meshgrid(x,x);

                nx = 2^(iter+1) + 1;
                dx = (domain(2)-domain(1))/(nx-1); x = domain(1):dx:domain(2);
                [X2 Y2] = meshgrid(x,x);

                if pars.IsCellX == 0
                    %x0 = interp2(X1,Y1, x0, X2,Y2,'linear');
                    x0 = interp2(X1,Y1, x0, X2,Y2,'spline');
                    %x0 = interp2(X1,Y1, x0, X2,Y2,'cubic');
                elseif pars.IsCellX == 1
                    for  di = 1: pars.nVarX
                        x0{di} = interp2(X1,Y1, x0{di}, X2,Y2,'cubic');
                    end
                end
        end

        % give boundary vaules to the initial point of level p+1
        x0 = feval(pars.boundary_value, x0, pars);
    end


    %     %----------------------------------------------------------------------
    %     % print optimal solution
    %     [nx, ny] = size(x_mg);
    %     dx = 1/(nx-1);
    %     x = 0:dx:1;
    %
    %     [X Y] = meshgrid(x,x);
    %     perf_fig = figure(iter+3);
    %     mesh(X, Y, x_mg)
    %     title(strcat('MG Solution', ' level = ' , num2str(pars.lev)) )
    %     print(perf_fig , '-depsc2', strcat('.\results\',pars.eval_fun,'_No_',...
    %         num2str(pars.model_index),'_D_',pars.direct_highlevel,...
    %         '_SM_',num2str(pars.iter_smooth),'_mg_lev',num2str(pars.lev_finest),'.eps'));
    %     %----------------------------------------------------------------------

end


output.pars = pars;

