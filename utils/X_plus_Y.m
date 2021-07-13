function z = X_plus_Y(x,y,alpha,pars)

% return z = x + alpha*y; 
% if there is boundary condition, return the same size of x 
if ~isequal(size(x),size(y)); error('the size of x and y are different!'); end;

%if pars.HaveBoundary == 1, the boundary of z will be set to zero

if pars.IsCellX == 0
    if pars.HaveBoundary == 0 % no boundary
       if alpha ~= 1 && alpha ~= -1; z = x + alpha*y; elseif alpha == 1; z = x + y; elseif alpha == -1; z = x - y; end
    elseif pars.HaveBoundary == 1 % dirichlit boundary
        if pars.Xinfo == 0;
            z = zeros(size(x)); 
            if alpha ~= 1 && alpha ~= -1;
                z(2:end-1) = x(2:end-1) + alpha*y(2:end-1);
            elseif alpha == 1;
                z(2:end-1) = x(2:end-1) + y(2:end-1);
            elseif alpha == -1;
                z(2:end-1) = x(2:end-1) - y(2:end-1);
            end
        elseif pars.Xinfo == 1;
            z = zeros(size(x)); 
            if alpha ~= 1 && alpha ~= -1;
                z(2:end-1,2:end-1) = x(2:end-1,2:end-1) + alpha*y(2:end-1,2:end-1);
            elseif alpha == 1;
                z(2:end-1,2:end-1) = x(2:end-1,2:end-1) + y(2:end-1,2:end-1);
            elseif alpha == -1;
                z(2:end-1,2:end-1) = x(2:end-1,2:end-1) - y(2:end-1,2:end-1);
            end            
        end
    end
elseif pars.IsCellX == 1
    z = cell(pars.nVarX,1);
    for di = 1:pars.nVarX
        if pars.HaveBoundary(di) == 0 % no boundary
            if alpha ~= 1 && alpha ~= -1; z{di} = x{di} + alpha*y{di}; 
            elseif alpha == 1; z{di} = x{di} + y{di}; 
            elseif alpha == -1; z{di} = x{di} - y{di}; end
        elseif pars.HaveBoundary(di) == 1  % dirichlit boundary
            if pars.Xinfo(di) == 0;
                z{di} = zeros(size(x{di})); 
                if alpha ~= 1 && alpha ~= -1;
                    z{di}(2:end-1) = x{di}(2:end-1) + alpha*y{di}(2:end-1);
                elseif alpha == 1;
                    z{di}(2:end-1) = x{di}(2:end-1) + y{di}(2:end-1);
                elseif alpha == -1;
                    z{di}(2:end-1) = x{di}(2:end-1) - y{di}(2:end-1);
                end
            elseif pars.Xinfo(di) == 1;
                z{di} = zeros(size(x{di}));
                if alpha ~= 1 && alpha ~= -1;
                    z{di}(2:end-1,2:end-1) = x{di}(2:end-1,2:end-1) + alpha*y{di}(2:end-1,2:end-1);
                elseif alpha == 1;
                    z{di}(2:end-1,2:end-1) = x{di}(2:end-1,2:end-1) + y{di}(2:end-1,2:end-1);
                elseif alpha == -1;
                    z{di}(2:end-1,2:end-1) = x{di}(2:end-1,2:end-1) - y{di}(2:end-1,2:end-1);
                end
            end
        end
    end
    
end


% 
% 
% % return z = x + alpha*y; 
% % if there is boundary condition, return the same size of x 
% if ~isequal(size(x),size(y)); error('the size of x and y are different!'); end;
% 
% %if pars.HaveBoundary == 1, the boundary of z will be set to the same as x
% z = x; 
% if pars.IsCellX == 0
%     if pars.HaveBoundary == 0 % no boundary
%        if alpha ~= 1 && alpha ~= -1; z = x + alpha*y; elseif alpha == 1; z = x + y; elseif alpha == -1; z = x - y; end
%     elseif pars.HaveBoundary == 1 % dirichlit boundary
%         if pars.Xinfo == 0;
%             if alpha ~= 1 && alpha ~= -1;
%                 z(2:end-1) = x(2:end-1) + alpha*y(2:end-1);
%             elseif alpha == 1;
%                 z(2:end-1) = x(2:end-1) + y(2:end-1);
%             elseif alpha == -1;
%                 z(2:end-1) = x(2:end-1) - y(2:end-1);
%             end
%         elseif pars.Xinfo == 1;
%             if alpha ~= 1 && alpha ~= -1;
%                 z(2:end-1,2:end-1) = x(2:end-1,2:end-1) + alpha*y(2:end-1,2:end-1);
%             elseif alpha == 1;
%                 z(2:end-1,2:end-1) = x(2:end-1,2:end-1) + y(2:end-1,2:end-1);
%             elseif alpha == -1;
%                 z(2:end-1,2:end-1) = x(2:end-1,2:end-1) - y(2:end-1,2:end-1);
%             end            
%         end
%     end
% elseif pars.IsCellX == 1
%     for di = 1:pars.nVarX
%         if pars.HaveBoundary(di) == 0 % no boundary
%             if alpha ~= 1 && alpha ~= -1; z{di} = x{di} + alpha*y{di}; 
%             elseif alpha == 1; z{di} = x{di} + y{di}; 
%             elseif alpha == -1; z{di} = x{di} - y{di}; end
%         elseif pars.HaveBoundary(di) == 1  % dirichlit boundary
%             if pars.Xinfo(di) == 0;
%                 if alpha ~= 1 && alpha ~= -1;
%                     z{di}(2:end-1) = x{di}(2:end-1) + alpha*y{di}(2:end-1);
%                 elseif alpha == 1;
%                     z{di}(2:end-1) = x{di}(2:end-1) + y{di}(2:end-1);
%                 elseif alpha == -1;
%                     z{di}(2:end-1) = x{di}(2:end-1) - y{di}(2:end-1);
%                 end
%             elseif pars.Xinfo(di) == 1;
%                 if alpha ~= 1 && alpha ~= -1;
%                     z{di}(2:end-1,2:end-1) = x{di}(2:end-1,2:end-1) + alpha*y{di}(2:end-1,2:end-1);
%                 elseif alpha == 1;
%                     z{di}(2:end-1,2:end-1) = x{di}(2:end-1,2:end-1) + y{di}(2:end-1,2:end-1);
%                 elseif alpha == -1;
%                     z{di}(2:end-1,2:end-1) = x{di}(2:end-1,2:end-1) - y{di}(2:end-1,2:end-1);
%                 end
%             end
%         end
%     end
%     
% end