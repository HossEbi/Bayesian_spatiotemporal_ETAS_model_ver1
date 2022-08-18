% Calculation of integral of rate of events with M=m in the interval [tstart,tend]
%Supplementary Figure S4 in Ebrahimian & Jalayer (2017) and also can be
% seen in Fig. 4 in Darzi et al. (2022)

% Mi & Ti  : the vectors of events

function N = calculate_N (Mi, Ti, ri, m, tstart, tend, Ml, K, theta, sum_intLambda, dA)%, xj, yj

beta  = theta(1);
c     = theta(2);
p     = theta(3);
d     = theta(4);
q     = theta(5);

Krt = (q-1)/pi*d^(2*(q-1))*(p-1)*c^(p-1);
% fun = @(x,y) 1./(( (x-xj)^2 + (y-yj)^2 + d^2 )^q);
% Kr = 1/integral2(fun,xj,inf,yj,inf);
% Krt = Kr*(p-1)*c^(p-1);

index = find(Ti < tend);    

if ~isempty(sum_intLambda) 
    
    if p == 1    
        Io = log((tend-Ti(index)+c)./(Ti(end)-Ti(index)+c));     
    else
        Io = ((tend-Ti(index)+c).^(1-p)-(Ti(end)-Ti(index)+c).^(1-p))/(1-p);
    end
    
    intLambda_end = exp(-beta*(m-Ml))*sum((K*exp(beta*(Mi(index)-Ml)).*Io).*dA.*(Krt./(ri(index).^2+d^2).^q));
    
    N = sum(exp(-beta*(m-Ml))*sum_intLambda) + intLambda_end;
    
else
    
    if p == 1
        Io = log((tend-Ti(index)+c)./(tstart-Ti(index)+c));     
    else
        Io = ((tend-Ti(index)+c).^(1-p)-(tstart-Ti(index)+c).^(1-p))/(1-p);
    end

    N = exp(-beta*(m-Ml))*sum((K*exp(beta*(Mi(index)-Ml)).*Io).*dA.*(Krt./(ri(index).^2+d^2).^q));
    
end

end