function TestBatchMGDriver

%clc;
%clear all; 
close all;

% change path
homeMGOPT = '/home/wenzw/work/MultigridIPM/code/MGLS-release';
cd(homeMGOPT);
addpath(strcat(homeMGOPT,'/Model'));
addpath(strcat(homeMGOPT,'/MG-OPT'));
addpath(strcat(homeMGOPT,'/MG-Lineq'));
addpath(strcat(homeMGOPT,'/utils'));

%--------------------------------------------------------------------------

% select a single test problem
Problist = {'minsurf1'}; % minimal surface
Problist = {'ell_nonlinear3'}; % nonlinear PDE
% Problist = {'nonconvex_ls'}; % the noncovex problem in Gratton, Sartenaer, Toint

%--------------------------------------------------------------------------

% choose a method
fullMG = 0;     % MLS-LBFGS
fullMG = 1;     % FMLS-LBFGS
% fullMG = 2;     % both MLS-LBFGS and FMLS-LBFGS
% fullMG = 3;     % MR-LBFGS
% fullMG = 4;     % FMLS-LMG
% fullMG = 5;     % FMLS-CG
% fullMG = 6;     % NT-FACT
% fullMG = 7;     % L-BFGS
fullMG = 14;    %  MLS-LBFGS, FMLS-LBFGS, MR-LBFGS
% fullMG = 15;    %  MLS-LBFGS, FMLS-LBFGS, MR-LBFGS, L-BFGS
% fullMG = 20;    %  MLS-LBFGS, FMLS-LBFGS, MR-LBFGS, FMLS-CG, L-BFGS, NT-FACT, FMLS-LMG 


% specify parameters
% mulitlevel information
finestlev = [8:10];  % finest level, could be a vector
% finestlev = [4:5];
finestlev = [10];
% finestlev = [8];

coarsestlev = 3;  % coarsest level, fixed scalar
pars.fine_lev_for_print = max(finestlev); % This parameter is only for print tables, appending empty lines

% print level
pars.display = 0;
% pars.display = 1;
% pars.display = 2;
% pars.display = 3;

% print figure or not
printfig = 0;
% printfig = 1;

% internal multigrid solver
pars.mg_opt = 'MGLS_Unc_Comp_L';

% tolerance of the norm of the graident
pars.eps_cvg = 1e-5;    %  decrease this value if there are too many direct search steps
                        % on the finest level or use "norm_g_low >= eps_cvg_low" in MGLS_Unc_Comp_L.m
pars.kappa_x = 0.1;
% pars.f_rel_tol = 1e-26;
% pars.tolx = 1e-12;

% % max iterations
% if strcmp(pars.enable_mg, 'yes') == 1
%     pars.maxiter = 100;
%     pars.maxiter_coarse = 3;
% elseif strcmp(pars.enable_mg, 'no') == 1
%     pars.maxiter = 200;
%     pars.maxiter_coarse = 1000;
% end

% pars.maxiter_coarse = 5;

% choose methods for the direct search
pars.direct_search = 'newton-fact';
% pars.direct_search = 'newton-pcg';
% pars.direct_search = 'newton-LMG';
pars.direct_search = 'lbfgs';

% number of storation for the lbfgs
pars.lbfgs_m = 5;

% inherit the information from last minimization sequence or not
pars.lbfgs_inherit_coarse_info = 0;
pars.lbfgs_inherit_coarse_info = 1;

pars.lbfgs_H_or_B = 'inverseB';
% pars.lbfgs_H_or_B = 'hessianB';

% initial point of the coarse level
pars.initial_x_coarse =  'old';
% pars.initial_x_coarse =  'restriction';

% number of smoothing
pars.iter_smooth = 0;
pars.iter_smooth = 1;

% pars.prolongation ='prolongation_regular';
% pars.restriction ='restriction_regular';

pars.int_op_2D = 'prolongate_bilinear2D';
pars.restrict_op_2D = 'restriction_bilinear2D';

pars.FMLS_int_op = 'bilinear';
pars.FMLS_int_op = 'cubic';

nProb = length(Problist);
for dprob = 1:nProb;
    %fprintf('\n\nTest problem: %s\n', Problist{di});
    % load problem
    pars.prob = Problist{dprob};
    pars = GetProb(Problist{dprob}, pars);

    % choose a method
    switch fullMG
        case {1, 2, 14, 15, 20}
            SolverVar = 'FMLS-LBFGS';
            for dfine = finestlev
                pars.enable_mg = 'yes';
                lev = dfine;
                pars.lev_coarest = coarsestlev;
                pars.lev_finest = lev;

                % choose methods for the direct search
                %pars.direct_search = 'newton-fact';
                % pars.direct_search = 'newton-pcg';
                % pars.direct_search = 'newton-LMG';
                pars.direct_search = 'lbfgs';
                
                %profile on;
                lev = pars.lev_coarest;
                nx = 2^lev + 1;
                if isscalar(pars.HaveBoundary)
                    u0 = zeros(nx,nx);
                    %u0 = 0.1*ones(nx,nx);
                else
                    u0{1,1} = zeros(nx,nx); u0{2,1} = zeros(nx,nx);
                    %u0{1,1} = rand(nx,nx); u0{2,1} = 10*randn(nx,nx);
                    
                    %dx = 1/(nx-1);  x = 0:dx:1; [X Y] = meshgrid(x,x); uu0 = sin(6*pi*X).*sin(2*pi*Y);
                    %u0{1,1} = uu0; u0{2,1} = zeros(nx,nx);
                end
                u0 = feval(pars.boundary_value, u0, pars);
                solve_t = cputime;
                [u_mg, output, stat] = Full_MGLS_Unc(u0, pars);
                output.cpu_time = cputime - solve_t;
                %profile viewer
                fprintf('\n--------------------------------------------------\n');
                fprintf('\nSolver: %s, Total Running Time: %f\n', SolverVar, output.cpu_time);
                OutputResults();
            end % end fineset level
    end % end switch

    switch fullMG

        % call multigrid method
        case {0, 2,14, 15, 20}
            SolverVar = 'MLS-LBFGS';
            for dfine = finestlev
                pars.enable_mg = 'yes';
                lev = dfine;
                pars.lev_coarest = coarsestlev;
                pars.lev_finest = lev;

                % choose methods for the direct search
                pars.direct_search = 'lbfgs';

                pars.iter_smooth = 1;

                %profile on;
                lev = pars.lev_finest;         nx = 2^lev + 1;
                if isscalar(pars.HaveBoundary)
                    u0 = zeros(nx,nx);
                    %u0 = 0.1*ones(nx,nx);
                else
                    u0{1,1} = zeros(nx,nx); u0{2,1} = zeros(nx,nx);
                    %u0{1,1} = 100*ones(nx,nx); u0{2,1} = u0{1,1};
                    
                    %dx = 1/(nx-1);  x = 0:dx:1; [X Y] = meshgrid(x,x); uu0 = sin(6*pi*X).*sin(2*pi*Y);
                    %u0{1,1} = uu0; u0{2,1} = zeros(nx,nx);

                end
                u0 = feval(pars.boundary_value, u0, pars);

                solve_t = cputime;
                [u_mg, output, stat] = MGLS_Unc(u0, pars);
                output.cpu_time = cputime - solve_t;
                fprintf('\n--------------------------------------------------\n');                
                fprintf('\nSolver: %s, Total Running Time: %f\n', SolverVar, output.cpu_time);

                %profile viewer
                OutputResults();

            end % end fineset level
    end % end switch

    switch fullMG
        % mesh-refinement
        case {3, 14, 15, 20}
            SolverVar = 'MR-LBFGS';
            for dfine = finestlev(end)

                pars.enable_mg = 'no';
                lev = dfine;
                pars.lev_coarest = coarsestlev;
                pars.lev_finest = lev;

                % choose methods for the direct search
                pars.direct_search = 'lbfgs';

                %profile on;
                lev = pars.lev_coarest;
                nx = 2^lev + 1;
                if isscalar(pars.HaveBoundary)
                    u0 = zeros(nx,nx);
                    %u0 = 0.1*ones(nx,nx);
                else
                    u0{1,1} = zeros(nx,nx); u0{2,1} = zeros(nx,nx);
                    %u0{1,1} = rand(nx,nx); u0{2,1} = 10*randn(nx,nx);
                end
                u0 = feval(pars.boundary_value, u0, pars);
                solve_t = cputime;
                [u_mg, output, stat] = Full_MGLS_Unc(u0, pars);
                output.cpu_time = cputime - solve_t;
                fprintf('\n--------------------------------------------------\n');                
                fprintf('\nSolver: %s, Total Running Time: %f\n', SolverVar, output.cpu_time);

                %profile viewer
                OutputResults();

            end % end fineset level
    end % end switch

    switch fullMG
        % multigrid, Newton's method using linear multigrid method
        case {4, 20}
            SolverVar = 'FMLS-LMG';
            for dfine = finestlev

                pars.enable_mg = 'yes';
                lev = dfine;
                pars.lev_coarest = coarsestlev;
                pars.lev_finest = lev;

                % choose methods for the direct search
                pars.direct_search = 'newton-LMG';

                %profile on;
                lev = pars.lev_coarest;
                nx = 2^lev + 1;
                if isscalar(pars.HaveBoundary)
                    u0 = zeros(nx,nx);
                else
                    u0{1,1} = zeros(nx,nx); u0{2,1} = zeros(nx,nx);
                    %u0{1,1} = rand(nx,nx); u0{2,1} = 10*randn(nx,nx);
                end
                u0 = feval(pars.boundary_value, u0, pars);
                solve_t = cputime;
                [u_mg, output, stat] = Full_MGLS_Unc(u0, pars);
                output.cpu_time = cputime - solve_t;
                fprintf('\n--------------------------------------------------\n');                
                fprintf('\nSolver: %s, Total Running Time: %f\n', SolverVar, output.cpu_time);

                %profile viewer
                OutputResults();

            end % end fineset level
    end % end switch

    switch fullMG
        % multigrid, Newton's method using PCG
        case {5, 20}
            SolverVar = 'FMLS-CG';
            for dfine = finestlev
                pars.enable_mg = 'yes';
                lev = dfine;
                pars.lev_coarest = coarsestlev;
                pars.lev_finest = lev;

                % choose methods for the direct search
                pars.direct_search = 'newton-pcg';

                %profile on;
                lev = pars.lev_coarest;
                nx = 2^lev + 1;
                if isscalar(pars.HaveBoundary)
                    u0 = zeros(nx,nx);
                else
                    u0{1,1} = zeros(nx,nx); u0{2,1} = zeros(nx,nx);
                    %u0{1,1} = rand(nx,nx); u0{2,1} = 10*randn(nx,nx);
                end
                u0 = feval(pars.boundary_value, u0, pars);
                solve_t = cputime;
                [u_mg, output, stat] = Full_MGLS_Unc(u0, pars);
                output.cpu_time = cputime - solve_t;
                fprintf('\n--------------------------------------------------\n');                
                fprintf('\nSolver: %s, Total Running Time: %f\n', SolverVar, output.cpu_time);

                %profile viewer
                OutputResults();

            end % end fineset level
    end % end switch

    switch fullMG
        % Newton's method using Factorization
        case {6, 20}
            SolverVar = 'NT-FACT';
            for dfine = finestlev
                pars.enable_mg = 'no';
                lev = dfine;
                pars.lev_coarest = lev;
                pars.lev_finest = lev;

                % choose methods for the direct search
                pars.direct_search = 'newton-fact';
                %profile on;
                lev = pars.lev_coarest;
                nx = 2^lev + 1;
                if isscalar(pars.HaveBoundary)
                    u0 = zeros(nx,nx);
                    %u0 = 0.1*ones(nx,nx);
                else
                    u0{1,1} = zeros(nx,nx); u0{2,1} = zeros(nx,nx);
                    %u0{1,1} = rand(nx,nx); 
                    %u0{2,1} = 0.1*ones(nx,nx);
                end
                u0 = feval(pars.boundary_value, u0, pars);
                solve_t = cputime;
                [u_mg, output, stat] = Full_MGLS_Unc(u0, pars);
                output.cpu_time = cputime - solve_t;
                fprintf('\n--------------------------------------------------\n');                
                fprintf('\nSolver: %s, Total Running Time: %f\n', SolverVar, output.cpu_time);

                %profile viewer
                OutputResults();

            end % end fineset level
    end % end switch

    switch fullMG
        % multigrid, Newton's method using linear multigrid method
        case {7, 15, 20}
            SolverVar = 'L-BFGS';
            for dfine = finestlev

                pars.enable_mg = 'no';
                lev = dfine;
                pars.lev_coarest = lev;
                pars.lev_finest = lev;

                % choose methods for the direct search
                pars.direct_search = 'lbfgs';

                %profile on;
                lev = pars.lev_coarest;
                nx = 2^lev + 1;
                if isscalar(pars.HaveBoundary)
                    u0 = zeros(nx,nx);
                else
                    u0{1,1} = zeros(nx,nx); u0{2,1} = zeros(nx,nx);
                    %u0{1,1} = rand(nx,nx); u0{2,1} = 10*randn(nx,nx);
                    
%                     dx = 1/(nx-1);  x = 0:dx:1; [X Y] = meshgrid(x,x);
%                     u0{1,1} = sin(6*pi*X).*sin(2*pi*Y);
                    
                end
                u0 = feval(pars.boundary_value, u0, pars);
                solve_t = cputime;
                [u_mg, output, stat] = Full_MGLS_Unc(u0, pars);
                output.cpu_time = cputime - solve_t;
                fprintf('\n--------------------------------------------------\n');                
                fprintf('\nSolver: %s, Total Running Time: %f\n', SolverVar, output.cpu_time);

                %profile viewer
                OutputResults();

            end % end fineset level
    end % end switch

end % end prob




% print and plot results
    function OutputResults()

        pars = output.pars;
        basename = strcat(pars.prob, SolverVar,  '_dir_',pars.direct_search,...
            '_MG_',pars.enable_mg, '_SM_',num2str(pars.iter_smooth), ...
            '_up_',num2str(pars.lev_finest), '_low_',num2str(pars.lev_coarest));

        filename = strcat('./results/',basename,'.txt');
        print_mg_output_spec(filename, pars, output, stat, SolverVar);

        %return;
        
        % perf_fig = figure(1);
        % mesh(X, Y, u_pde)
        % title('True Solution');
        % print(perf_fig , '-depsc', strcat('.\results\',pars.eval_fun,'_No_',...
        %         num2str(pars.model_index),'_exact_lev',num2str(pars.lev_finest),'.eps'));
        %
        %
        domain = pars.domain;
        nx = 2^pars.lev_finest + 1;     ny = nx;
        dx = (domain(2)-domain(1))/(nx-1);    dy = (domain(4)-domain(3))/(ny-1);
        area = 0.5*dx*dy;    x = domain(1):dx:domain(2);    [X Y] = meshgrid(x,x);

        ifig = 1;
        if pars.IsCellX == 0
            ifig = ifig + 1; perf_fig = figure(ifig); mesh(X, Y, u_mg);
            title(strcat('MG Solution', ' level = ' , num2str(pars.lev_finest)) );
            %print(perf_fig , '-depsc', strcat('.\results\',basename,'.eps'));
        elseif pars.IsCellX == 1;
            for di = 1: pars.nVarX
                ifig = ifig + 1; perf_fig = figure(ifig); mesh(X, Y, u_mg{di});
                title(strcat('MG Solution', ' level = ' , num2str(pars.lev_finest)) );
%                 if printfig == 1
%                     print(perf_fig , '-depsc2', strcat('.\results\', basename,'_solu_x',num2str(di),'.eps'));
%                 end
            end
        end


        %perf_fig = figure(5);    mesh(X, Y, u_pde);
        ifig = ifig + 1; perf_fig = figure(ifig);
        yp = log10(stat.vec_norm_g);
        nit = length(yp); idx = 1:nit; idmgr = logical(stat.x_type); idmg = ~logical(stat.x_type);
        plot(idx(idmgr), yp(idmgr), 'or', 'LineWidth',2, 'MarkerSize',7);
        hold on;
        plot(idx(idmg), yp(idmg), 'x', 'LineWidth',2, 'MarkerSize',4);
        hold off;
        legend('recursive', 'direct', 0);
        xlabel('index'); ylabel('log10(norm(g))');
        title('norm of the gradient');
        if ~isempty(findstr(SolverVar, 'FMLS') )  %fullMG ~= 0
            hold on;
            ndi = 20;
            xmin = 0; xmax = nit+10; ymin = min(yp); ymax = max(yp);
            ac1 = 1.1; ac2 = 0.95;
            if ymin < 0; ymin = ac1*ymin; else ymin = ac2*ymin; end
            if ymax < 0; ymax = ac2*ymax; else ymax = ac1*ymax; end
            for di = pars.lev_coarest+1:pars.lev_finest-1
                plot((stat.iter_gf(di) + 0.5) *ones(ndi,1), ymin:(ymax-ymin)/(ndi-1):ymax, '-.');
            end
            hold off;    %axis([xmin xmax ymin ymax]);
        end
        if printfig == 1
            print(perf_fig , '-depsc2', strcat('.\results\', basename,'_grad','.eps'));
        end

        ifig = ifig + 1; perf_fig = figure(ifig);
        yp = log10(abs(stat.vec_f));
        nit = length(yp); idx = 1:nit; idmgr = logical(stat.x_type); idmg = ~logical(stat.x_type);
        %plot(idx(idmgr), yp(idmgr), 'o', idx(idmg), yp(idmg), 'x', 'LineWidth',3,'MarkerSize',5);
        plot(idx(idmgr), yp(idmgr), 'or', 'LineWidth',2, 'MarkerSize',7);
        hold on;
        plot(idx(idmg), yp(idmg), 'x', 'LineWidth',2, 'MarkerSize',4);
        hold off;
        legend('recursive', 'direct', 0);
        xlabel('index'); ylabel('log10(abs(f))');
        title('function value f');
        if  ~isempty(findstr(SolverVar, 'FMLS') )  %fullMG == 1
            hold on;
            ndi = 20;
            xmin = 0; xmax = nit+10; ymin = min(yp); ymax = max(yp);
            ac1 = 1; ac2 = 1;
            if ymin < 0; ymin = ac1*ymin; else ymin = ac2*ymin; end
            if ymax < 0; ymax = ac2*ymax; else ymax = ac1*ymax; end
            for di = pars.lev_coarest+1:pars.lev_finest-1
                plot((stat.iter_gf(di) + 0.5) *ones(ndi,1), ymin:(ymax-ymin)/(ndi-1):ymax, '-.'); %,'LineWidth',2, 'MarkerSize',8
            end
            hold off;    %axis([xmin xmax ymin ymax]);
        end
        if printfig == 1
            print(perf_fig , '-depsc2', strcat('.\results\', basename,'_fun','.eps'));
        end

        ifig = ifig + 1; perf_fig = figure(ifig);
        plot([1:stat.iter_count], stat.iter_hist(1:stat.iter_count),'*', ...
            [1:stat.iter_count], stat.iter_hist(1:stat.iter_count),'-');

        axis([0 stat.iter_count+1 1 pars.lev_finest+1]);
        title('Iteration History');
        xlabel('Iteration');    ylabel('level');
        if printfig == 1
            print(perf_fig , '-depsc2', strcat('.\results\', basename,'_Ithist','.eps'));
        end

        %save(strcat('.\results\', basename));
        %save('minsurf_eps_1000');
        %fprintf('||u_mg - u_pde||: %3.2e, rel.err: %3.2e \n', norm(u_mg - u_pde, 'fro'), norm(u_mg - u_pde, 'fro')/norm(u_mg, 'fro'));


    end


end