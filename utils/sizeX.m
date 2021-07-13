function z = sizeX(x,pars)

if pars.IsCellX == 0
    if pars.HaveBoundary == 0 % no boundary
        [mx,nx] = size(x);   z = mx*nx;
    elseif pars.HaveBoundary == 1 % dirichlit boundary
        if pars.Xinfo == 0;
            z = length(x); if z < 2; error('length of x is %d < 2', z); else z = z -2; end;
        elseif pars.Xinfo == 1;
            [mx,nx] = size(x);   
            if mx < 2 || nx < 2; error('size of x is [%d, %d] < 2', mx,nx); else z = (mx-2)*(nx-2); end;
        end
    end
elseif pars.IsCellX == 1
    z = 0;
    for di = 1:pars.nVarX
        if pars.HaveBoundary(di) == 0 % no boundary
            [mx,nx] = size(x{di});   z = z + mx*nx;
        elseif pars.HaveBoundary(di) == 1  % dirichlit boundary
            if pars.Xinfo == 0;
                n = length(x{di}); if z < 2; error('length of x is %d < 2', z); else z = z+n -2; end;
            elseif pars.Xinfo == 1;
                [mx,nx] = size(x{di});
                if mx < 2 || nx < 2; error('size of x is [%d, %d] < 2', mx,nx); else z = z+(mx-2)*(nx-2); end;
            end
        end
    end
end