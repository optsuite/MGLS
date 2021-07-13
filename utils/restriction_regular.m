function u_l = restriction_regular(u_p, pars)

%--------------------------------------------------------------------------
% Prolongation operator for problem with Dirichlet boundary condition,
%   then, you may enforce nodes on the boundary satisfy the boundary
%   condtion.

% input:    
%   u_l , value on lower level
%   u_p , the boundary has been initialized 
%
% Nine points restriction
    
if pars.IsCellX == 0
       if pars.Xinfo == 0;
           error('the 1-dimensional prolongation has not implemented yet');
       elseif pars.Xinfo == 1;
            if strcmp(pars.var_type, 'dir') == 1
                u_l = feval(pars.restrict_op_2D, u_p,pars);
            elseif strcmp(pars.var_type, 'var') == 1
                u_l = restriction_simple2D(u_p,pars);
            end
       end
elseif pars.IsCellX == 1
    u_l = cell(pars.nVarX,1);
    for di = 1:pars.nVarX
       if pars.Xinfo(di) == 0;
           error('the 1-dimensional prolongation has not implemented yet');
       elseif pars.Xinfo(di) == 1;
            if strcmp(pars.var_type, 'dir') == 1
                u_l{di} = feval(pars.restrict_op_2D, u_p{di},pars);
            elseif strcmp(pars.var_type, 'var') == 1
                u_l{di} = restriction_simple2D(u_p{di},pars);
            end            
       end
    end
end

end

