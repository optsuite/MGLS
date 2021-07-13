function z = X_minus_Y(x,y,pars)

% return z = x + y; 
% if there is boundary condition, return the same size of x 
if ~isequal(size(x),size(y)); error('the size of x and y are different!'); end;

if pars.IsCellX == 0
    if pars.HaveBoundary == 0 % no boundary
       z = x - y; 
    elseif pars.HaveBoundary == 1 % dirichlit boundary
        if pars.Xinfo == 0;
            z = zeros(size(x)); z(2:end-1) = x(2:end-1) - y(2:end-1);
        elseif pars.Xinfo == 1;
            z = zeros(size(x)); z(2:end-1,2:end-1) = x(2:end-1,2:end-1) - y(2:end-1,2:end-1);
        end
    end
elseif pars.IsCellX == 1
    z = cell(pars.nVarX,1);
    for di = 1:pars.nVarX
        if pars.HaveBoundary(di) == 0 % no boundary
            z{di} = x{di} - y{di};
        elseif pars.HaveBoundary(di) == 1  % dirichlit boundary
            if pars.Xinfo(di) == 0;
                z{di} = zeros(size(x{di})); z{di}(2:end-1) = x{di}(2:end-1) - y{di}(2:end-1);
            elseif pars.Xinfo(di) == 1;
                z{di} = zeros(size(x{di})); z{di}(2:end-1,2:end-1) = x{di}(2:end-1,2:end-1) - y{di}(2:end-1,2:end-1);
            end
        end
    end
    
end