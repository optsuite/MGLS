function z = normX(x,mode,pars)

if pars.IsCellX == 0
    if pars.HaveBoundary == 0 % no boundary
        z = norm(x,mode);
    elseif pars.HaveBoundary == 1 % dirichlit boundary
        if pars.Xinfo == 0;
            z = norm(x(2:end-1),mode);
        elseif pars.Xinfo == 1;
            z = norm(x(2:end-1,2:end-1), mode);
        end
    end
    
%     [nx,ny] = size(x);     dx = 1/(nx-1);    dy = 1/(ny-1);     area = 0.5*dx*dy;
%     z = z*sqrt(area);

elseif pars.IsCellX == 1
    z = zeros(pars.nVarX,1);
    for di = 1:pars.nVarX
        if pars.HaveBoundary(di) == 0 % no boundary
            z(di) = norm(x{di}, mode);
        elseif pars.HaveBoundary(di) == 1  % dirichlit boundary
            if pars.Xinfo(di) == 0;
                z(di) = norm(x{di}(2:end-1),mode);
            elseif pars.Xinfo(di) == 1;
                z(di) = norm(x{di}(2:end-1,2:end-1),mode);
            end
        end
    end

    if isequal(mode,1)
        z= sum(z);
    elseif isequal(mode,inf)
        z = max(z);
    elseif isequal(mode,2) || isequal(mode,'fro')
        z =  sqrt(sum(z.^2));
    end
    
%     [nx,ny] = size(x{1});     dx = 1/(nx-1);    dy = 1/(ny-1);     area = 0.5*dx*dy;
%     z = z*sqrt(area);
end

