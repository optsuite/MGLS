function [z, sizex] = VecX(x,pars)

if pars.IsCellX == 0
    if pars.HaveBoundary == 0 % no boundary
        sizex = size(x);
        if pars.Xinfo == 1;
            z = x(:);
        end        
    elseif pars.HaveBoundary == 1 % dirichlit boundary
        [mX,nX] = size(x); sizex = [mX,nX];
        if pars.Xinfo == 1;
            z = reshape(x(2:end-1,2:end-1), (mX-2)*(nX-2),1);
        end
    end
elseif pars.IsCellX == 1
    z = []; sizex = cell(pars.nVarX,1);
    for di = 1:pars.nVarX
        if pars.HaveBoundary(di) == 0 % no boundary
            sizex{di} = size(x{di});
            if pars.Xinfo == 1;
                z = [z; x{di}(:)];
            end
        elseif pars.HaveBoundary(di) == 1  % dirichlit boundary
            [mX,nX] = size(x{di}); sizex{di} = [mX,nX];
            if pars.Xinfo == 1;
                z = [z; reshape(x{di}(2:end-1,2:end-1), (mX-2)*(nX-2),1)];
            end
        end
    end
end