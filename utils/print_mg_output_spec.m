function print_mg_output_spec(filename,pars, output, stat, solver)


output

lev = (pars.lev_coarest:pars.lev_finest)';
stat.iter = stat.iter(lev);
stat.nfe = stat.nfe(lev);
stat.nge = stat.nge(lev);
stat.nhe = stat.nhe(lev);
stat.n_cycle = stat.n_cycle(lev);
stat.norm_g = stat.norm_g(lev);
stat.cpu = stat.cpu(lev);


fprintf('\n\nModel Name: %s\n', pars.eval_fun);
fprintf('Solver: %s on level %d\n\n', solver, pars.lev_finest);
fprintf('lev \t nls \t nfe \t nge \t nhe \t nv \t ||g|| \t CPU\n');
fprintf('%d \t %d \t %d \t %d \t %d \t %d \t %3.1e \t %6.2f\n', [lev stat.iter stat.nfe stat.nge stat.nhe stat.n_cycle  stat.norm_g stat.cpu]');

fid = fopen(filename,'w');

fprintf(fid,'Model Name: %s\n\n', pars.eval_fun);

% print out table
fprintf(fid,'\\begin{tabular}[t]{|');

item = pars.lev_finest -pars.lev_coarest +1;
for iter = 1:item+1
   fprintf(fid,'c| ' ); 
end
fprintf(fid,'}\\hline \n');
fprintf(fid,' & \\multicolumn{%d}{|c|}{%s on level %d} \\\\ \\hline \n', item, solver, pars.lev_finest);
fprintf(fid,'\t h & \t');
for iter = pars.lev_coarest : pars.lev_finest-1
   fprintf(fid,'$%d$ &\t', iter ); 
end
fprintf(fid,'$%d$ \\\\ \\hline \n', pars.lev_finest );

printout_vec(fid,pars, 'nls', stat.iter);
printout_vec(fid,pars, 'nfe', stat.nfe);
printout_vec(fid,pars, 'nge', stat.nge);
%printout_vec(fid,pars, 'nhe', stat.nhe);
printout_vec(fid,pars, 'ncycle', stat.n_cycle);
printout_vec_double(fid,pars, '$||g^*||$', stat.norm_g);
printout_vec_double(fid,pars, 'CPU', stat.cpu);
%fprintf(fid,'$ \\|g_{\\mathrm{N}}^* \\|_2 $ & \\multicolumn{%d}{|c|}{%3.2e} \\\\ \\hline \n', item, output.norm_g);
%fprintf(fid,'CPU & \\multicolumn{%d}{|c|}{%3.2e} \\\\ \\hline \n', item, output.cpu_time);
fprintf(fid,'\\end{tabular}\n');


app = zeros(pars.fine_lev_for_print - pars.lev_finest,1);
stat.iter = [stat.iter; app];
stat.nfe = [stat.nfe; app];
stat.nge = [stat.nge; app];
stat.nhe = [stat.nhe; app];
stat.n_cycle = [stat.n_cycle; app];
stat.norm_g = [stat.norm_g; app];
stat.cpu = [stat.cpu; app];
lev = (pars.lev_coarest:pars.fine_lev_for_print)';

fprintf(fid,'\n\n');
sstat = [lev,  stat.nfe, stat.nge, stat.n_cycle, stat.norm_g,stat.cpu];
fprintf(fid,'\\begin{tabular}[t]{|c|c|c|c|c|c|}\\hline \n');
fprintf(fid,'\\multicolumn{6}{|c|}{%s on level %d} \\\\ \\hline \n', solver, pars.lev_finest);
fprintf(fid,'h & nfe & nge & nv & $||g^*||$ & CPU \\\\ \\hline \n');
fprintf(fid,'%d & %d & %d & %d & %3.1e & %6.2f \\\\ \\hline \n', sstat');
fprintf(fid,'\\end{tabular}\n');

fprintf(fid,'\n\nwithout n_cycle \n');
sstat = [lev,  stat.nfe, stat.nge, stat.norm_g,stat.cpu];
fprintf(fid,'\\begin{tabular}[t]{|c|c|c|c|c|}\\hline \n');
fprintf(fid,'\\multicolumn{5}{|c|}{%s on level %d} \\\\ \\hline \n', solver, pars.lev_finest);
fprintf(fid,'h & nfe & nge &  $||g^*||$ & CPU \\\\ \\hline \n');
fprintf(fid,'%d & %d & %d &  %3.1e & %6.2f \\\\ \\hline \n', sstat');
fprintf(fid,'\\end{tabular}\n');


fprintf(fid,'\n\n MLS \n\n');
sstat = [lev,  stat.nfe, stat.nge, stat.n_cycle];
fprintf(fid,'\\begin{tabular}[t]{|c|c|c|c|}\\hline \n');
fprintf(fid,'\\multicolumn{4}{|c|}{%s on level %d} \\\\ \\hline \n', solver, pars.lev_finest);
fprintf(fid,'h & nfe & nge &nv  \\\\ \\hline \n');
fprintf(fid,'%d & %d & %d & %d \\\\ \\hline \n', sstat');
fprintf(fid,'$ ||g^*||$ & \\multicolumn{3}{|c|}{%3.2e} \\\\ \\hline \n', stat.norm_g(pars.lev_finest-pars.lev_coarest +1));
fprintf(fid,'CPU & \\multicolumn{3}{|c|}{%6.2f} \\\\ \\hline \n', stat.cpu(pars.lev_finest-pars.lev_coarest +1));
fprintf(fid,'\\end{tabular}\n');


% fprintf(fid,'\n\n');
% sstat = [(pars.lev_coarest: pars.lev_finest)', stat.iter , stat.nfe, stat.nge, stat.n_cycle, stat.norm_g,stat.cpu];
% fprintf(fid,'\\begin{tabular}[t]{|c|c|c|c|c|c|c|}\\hline');
% fprintf(fid,'\\multicolumn{7}{|c|}{a} \\\\ \\hline \n');
% fprintf(fid,'h & nls & nfe & nge & nv & $||g^*||$ & CPU \\\\ \\hline \n');
% fprintf(fid,'%d & %d & %d & %d & %d & %3.1e & %3.1e \\\\ \\hline \n', sstat');
% fprintf(fid,'\\end{tabular}\n');

PrintStructTree(output, fid);
PrintStructTree(pars, fid);

fclose(fid);

% print out statistics 
function printout_vec(fid,pars, name, vec)
fprintf(fid,'\t %s & \t', name);
for iter = 1 : pars.lev_finest-pars.lev_coarest
    fprintf(fid,'%d & \t', vec(iter) ); 
end
fprintf(fid,'%d \\\\ \\hline \n', vec(iter+1) );

function printout_vec_double(fid,pars, name, vec)
fprintf(fid,'\t %s & \t', name);
for iter = 1 : pars.lev_finest-pars.lev_coarest
    fprintf(fid,'%3.1e & \t', vec(iter) ); 
end
fprintf(fid,'%3.1e \\\\ \\hline \n', vec(iter+1) );