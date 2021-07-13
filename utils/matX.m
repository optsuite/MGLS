function z = matX(x,sizex,pars)

if pars.IsCellX == 0
    if pars.HaveBoundary == 0 % no boundary
        z = reshape(x, sizex);
    elseif pars.HaveBoundary == 1 % dirichlit boundary
        if pars.Xinfo == 0;
            z = zeros(sizex); z(2:end-1) = x(2:end-1);
        elseif pars.Xinfo == 1;            
            z = zeros(sizex); z(2:end-1,2:end-1) = reshape(x, sizex-2);
        end
    end
elseif pars.IsCellX == 1
    z = cell(pars.nVarX,1); idx = 1;
    for di = 1:pars.nVarX
        if pars.HaveBoundary(di) == 0 % no boundary
            t = prod(sizex{di});
            z{di} = reshape(x(idx:idx+t-1), sizex{di} );
            idx = idx+t;
        elseif pars.HaveBoundary(di) == 1  % dirichlit boundary
            if pars.Xinfo == 0;
                t = prod(sizex{di}(1))-2;
                z{di} = zeros(sizex{di}); z{di}(2:end-1) = x(idx:idx+t-1);
                idx = idx+t;
            elseif pars.Xinfo == 1;
                %z{di} = zeros(sizex{di}); z{di}(2:end-1,2:end-1) = reshape(x{di}, sizex{di}-2);
                z{di} = zeros(sizex{di}); t = prod(sizex{di}-2);
                z{di}(2:end-1,2:end-1) = reshape(x(idx:idx+t-1), sizex{di}-2 );
                idx = idx+t;

            end
        end
    end
end