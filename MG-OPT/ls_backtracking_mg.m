function [stp, f, options, work] = ls_backtracking_mg(stp,f, g, forigin, deriv_cd_g, options , work)
%--------------------------------------------------------------------------
%
% Line search with backtracking
% Given g' d < 0 and ftol < 0.5, find xnew = x + a d, such that
%   f(xnew) < f(x) + a * ftol * g' d, using a backtracking line search
% and f - forigin >= (1 - ftol)* g0' ( x + a*d - x0  )
% 
% Reference:
%   J.E. Dennis, JR, Robert B. Schnabel
%   Numerical Methods for Unconstrained Optimization
%     and Nonlinear Equations
%   Page 325-327
%
%
% Author: Zaiwen Wen
%   2007-10-15
%--------------------------------------------------------------------------

% Inverse Communication

zero=0.0;    

if  work.task == 1

    % default options
    if ~isfield(options,'maxiter')
        options.maxiter = 20;
    end

    if ~isfield(options,'display')
        options.display = 'no';
%         options.display = 'iter';
    end

    % options for wolfe condition
    if ~isfield(options,'ftol')
        options.ftol = 1e-3;
    end

    if ~isfield(options,'xtol')
        options.xtol = 1e-30;
    end

    if ~isfield(options,'stpmin')
        options.stpmin = 1e-20;
        %         options.stpmin = 1e-10;
    end

    if ~isfield(options,'stpmax')
        options.stpmax = 1e5;
        %         options.stpmax = 2;
    end

    % c        Check the input arguments for errors.
    if (stp < options.stpmin)
        work.task = -5;
        work.msg = 'ERROR: STP .LT. STPMIN';
    end

    if (stp > options.stpmax)
        work.task = -5;
        work.msg = 'ERROR: STP .GT. STPMAX';
    end

    if (g > zero)
        work.task = -5;
        work.msg = 'ERROR: INITIAL G .GE. ZERO';
    end

    if (options.ftol < zero)
        work.task = -5;
        work.msg =  'ERROR: FTOL .LT. ZERO';
    end

    if (options.xtol < zero)
        work.task = -5;
        work.msg =  'ERROR: XTOL .LT. ZERO';
    end

    % c        Exit if there are errors on input.
    if work.task == -5
        error(work.msg);
        return;
    end

    % c        Initialize local variables.

    % c        The variables stx, fx, gx contain the values of the step,
    % c        function, and derivative at the best step.
    % c        The variables sty, fy, gy contain the value of the step,
    % c        function, and derivative at sty.
    % c        The variables stp, f, g contain the values of the step,
    % c        function, and derivative at stp.

    work.initslope = g;
    work.gtest     = options.ftol * work.initslope;
    work.finit     = f;
   
    work.bestfx    = f;
    work.beststp   = stp;

    work.iter = 0;
    work.task = 2;   
    work.msg = 'FG';

    return;

end

%------------------------------------------------------
% main loop
%  if  work.task == 2, the code come back with new f and g
if  work.task == 2 &&  work.bestfx > f
    work.bestfx    = f;
    work.beststp   = stp;
end

% exceed max iterations
if (work.iter >= options.maxiter)
    % currently, just set work.task = 0
    work.task = 0;
    work.msg =  'EXCEED MAX ITERATIONS';

    % set the stp to best in the history
    stp = work.beststp;
    f   = work.bestfx;
    return;
end

work.iter = work.iter + 1;

gtol = 1 - options.ftol;
% gtol = 0.2;

% ftest = work.finit + stp*work.gtest;
% if strcmp(options.display, 'iter')
%    fprintf('iter: %2d \t stp: %3.2e \t crit1: %3.2e \t crit2: %3.2e \n', ...
%             work.iter, stp, f - ftest, -f + forigin + (1 - options.ftol)*deriv_cd_g); 
% end
% 
% % c     Test for convergence.
% if f <= ftest && f - forigin - (1 - options.ftol)*deriv_cd_g > 0
%     work.task = 0;
%     work.msg =  'CONVERGENCE';
%     return;
% end

ftest = work.finit + stp*work.gtest;
if strcmp(options.display, 'iter')
   fprintf('iter: %2d \t stp: %3.2e \t crit1: %3.2e \t crit2: %3.2e \n', ...
            work.iter, stp, f - ftest, -f + forigin + gtol*deriv_cd_g); 
end

% c     Test for convergence.
if f <= ftest && f - forigin - gtol*deriv_cd_g > 0
    work.task = 0;
    work.msg =  'CONVERGENCE';
    return;
end

beta = 0.1;
stp = beta*stp;
% c     Obtain another function and derivative.
work.task = 2;     work.msg = 'FG';
return;


% % c     Test for warnings.
% if stp >= options.stpmax
%     work.task = -1;
%     work.msg =  'WARNING: STP = STPMAX';
% end
% 
% if stp <= options.stpmin
%     work.task = -1;
%     work.msg =  'WARNING: STP = STPMIN';
% end
% 
% % if (task(1:4) .eq. 'WARN' .or. task(1:4) .eq. 'CONV') go to 10
% if (work.task == -1 )
%     stp = work.beststp;
%     f   = work.bestfx;
%     return;
% end
% 
% % the first backtrack, quadratic fit
% if work.iter == 1
%     stptemp = -work.initslope/ ( 2* (f - work.finit - work.initslope));
% % all subsequent backtracks, cubic fit
% else
%     coeff = 1/(stp - work.stpprev) * ...
%             [ 1/stp^2,              -1/(work.stpprev.^2);   
%              -work.stpprev/ stp^2,   stp/(work.stpprev.^2)] * ...
%              [ f - work.finit - stp* work.initslope;
%                work.fprev - work.finit - work.stpprev * work.initslope];
%     disc = coeff(2)^2 - 3* coeff(1)* work.initslope;
%     
%     if disc < 0 
%         % erorr('disc < 0');
%         % disc = 0;
%         stp = 0.1*stp;
%         return;
%     end
%     
%     if coeff(1) == 0
%         stptemp = - work.initslope/(2*coeff(2));
%     else
%         stptemp = ( -coeff(2) + sqrt(disc)) / (3*coeff(1));
%     end
%              
%     if stptemp > 0.5*stp
%         stptemp = 0.5*stp;
%     end
% end
% 
% if strcmp(options.display, 'iter')
%    fprintf('iter: %2d \t stp: %3.2e \t stptemp: %3.2e \n', work.iter, stp, stptemp); 
% end
% 
% work.stpprev    = stp;
% work.fprev      = f;
% 
% if stptemp <= 0.1* stp
%     stp = 0.1*stp;
% else
%     stp = stptemp;
% end
% 
% 
% % c     Obtain another function and derivative.
% work.task = 2;     work.msg = 'FG';