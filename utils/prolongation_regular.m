function u_p = prolongation_regular(u_l,  pars)

%--------------------------------------------------------------------------
% Prolongation operator for problem with Dirichlet boundary condition,
%   then, you may enforce nodes on the boundary satisfy the boundary
%   condtion.

% input:    
%   u_l , value on lower level
%   u_p , the boundary has been initialized 

if pars.IsCellX == 0
       if pars.Xinfo == 0;
           error('the 1-dimensional prolongation has not implemented yet');
       elseif pars.Xinfo == 1;
            u_p = feval(pars.int_op_2D, u_l,pars);
       end
elseif pars.IsCellX == 1
    u_p = cell(pars.nVarX,1);
    for di = 1:pars.nVarX
       if pars.Xinfo(di) == 0;
           error('the 1-dimensional prolongation has not implemented yet');
       elseif pars.Xinfo(di) == 1;
            u_p{di} = feval(pars.int_op_2D,u_l{di},pars);
       end
    end
end

end
