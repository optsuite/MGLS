function z = ddotXY(x,y,alpha,pars)

% return z = x + alpha*y; 
% if there is boundary condition, return the same size of x 
if ~isequal(size(x),size(y)); error('the size of x and y are different!'); end;

%if pars.HaveBoundary == 1, the boundary of z will be set to zero

if pars.IsCellX == 0
    if pars.HaveBoundary == 0 % no boundary
        if pars.Xinfo == 0;
            z = dot(x, y);
        elseif pars.Xinfo == 1;
            z = sum(sum(x.*y));
        end
    elseif pars.HaveBoundary == 1 % dirichlit boundary
        if pars.Xinfo == 0;
            z = dot(x(2:end-1), y(2:end-1));
        elseif pars.Xinfo == 1;
            z = sum(sum(x(2:end-1,2:end-1).*y(2:end-1,2:end-1)));
        end
    end
elseif pars.IsCellX == 1
    z = 0;
    for di = 1:pars.nVarX
        if pars.HaveBoundary(di) == 0 % no boundary
            if pars.Xinfo == 0;
                z = z+dot(x{di}, y{di});
            elseif pars.Xinfo == 1;
                z = z+sum(sum(x{di}.*y{di}));
            end
        elseif pars.HaveBoundary(di) == 1  % dirichlit boundary
            if pars.Xinfo == 0;
                z = z + dot(x{di}(2:end-1), y{di}(2:end-1));
            elseif pars.Xinfo == 1;
                z = z + sum(sum(x{di}(2:end-1,2:end-1).*y{di}(2:end-1,2:end-1)));
            end
        end
    end
end
z = alpha*z;