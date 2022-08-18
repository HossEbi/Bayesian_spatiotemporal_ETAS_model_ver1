function [y] = interpola(X,Y,x)

if (x > min(X) && x < max(X))
    
    indupp = find(X > x, 1, 'first');
    Xmax = X(indupp);
    
    indlow = find(X < x, 1, 'last');
    Xmin = X(indlow);
    
    if (Xmax == Xmin)
        disp('Xmin is equal to Xmax')
        y = Y(indupp);
    else
        y = interp1([Xmin,Xmax],[Y(indlow),Y(indupp)],x);        
    end
    
elseif x <= min(X)
    
    disp('warning! x is lower than Xmin')
    y = interp1([X(1),X(2)],[Y(1),Y(2)],x,'linear','extrap');      
    
elseif x >= max(X)    
    
    disp('warning! x is higher than Xmax')
    y = interp1([X(end-1),X(end)],[Y(end-1),Y(end)],x,'linear','extrap');  
    
end

