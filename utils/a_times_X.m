function z = a_times_X(alpha,x,pars)

if pars.IsCellX == 0
    z = alpha*x;
elseif pars.IsCellX == 1
    z = x;
    for di = 1:pars.nVarX
        z{di} = alpha*x{di};
    end
end