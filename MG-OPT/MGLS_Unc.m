function [x, output, stat] = MGLS_Unc(x0, pars,  varargin)

% The dirver file for the multigrid method
%
%--------------------------------------------------------------------------
% Copyright (c) 2008.  Zaiwen Wen 
%--------------------------------------------------------------------------


if isempty(x0)
    error('please provide an initial solution!');
end

pars = MGLS_Unc_pars_set(x0, pars);
lev = pars.lev_finest;         
pars.lev = pars.lev_finest;

stat.nfe = zeros(lev,1);
stat.nge = zeros(lev,1);
stat.nhe = zeros(lev,1);
stat.norm_g = zeros(lev,1);
stat.cpu = zeros(lev,1);
stat.n_cycle = zeros(lev,1);
stat.iter = zeros(lev,1);
stat.iter_count = 0;
stat.iter_hist = zeros(10^3,1);

% if possible, build projection and prolongation matrix
if strcmp(pars.direct_search, 'newton-LMG') == 1
    pars = build_pro_restr_matrix( pars);
    pars.have_pro_restr = 1;
end

tic;
[x, x0, output, stat] = feval(pars.mg_opt, x0, [], pars, stat, varargin{:});
stat.cpu(pars.lev) = toc;
stat.norm_g(pars.lev) = output.norm_g;
stat.iter_gf(pars.lev) = length(stat.vec_norm_g);
output.pars = pars;
