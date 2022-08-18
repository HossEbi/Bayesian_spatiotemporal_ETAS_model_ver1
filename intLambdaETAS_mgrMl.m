% Function for calculating the integral of Lambda_ETAS(t,m>Ml) from tstart to tend
% Mi & Ti  : the vectors of events
% written by: Hossein Ebrahimian   
% Last update: 01/2017

function intLambda = intLambdaETAS_mgrMl (Mi, Ti, tstart, tend, Ml, K, theta)
%tstart = 0.0; tend = time(1);K = 1;

beta  = theta(1);
c     = theta(2);
p     = theta(3);

Kt = (p-1)*c^(p-1);

index = find(Ti < tend); 

if p == 1    
    Io = log((tend-Ti(index)+c)./(tstart-Ti(index)+c));     
else    
%     if (tend-Ti(index)+c)>=0 && (tstart-Ti(index)+c)>=0
        Io = ((tend-Ti(index)+c).^(1-p)-(tstart-Ti(index)+c).^(1-p))/(1-p);
%     elseif (tend-Ti(index)+c)>=0 && (tstart-Ti(index)+c)<=0
%         Io = ((tend-Ti(index)+c).^(1-p)-(-1).*abs(tstart-Ti(index)+c).^(1-p))/(1-p);                
%     elseif (tend-Ti(index)+c)<=0 && (tstart-Ti(index)+c)>=0
%         Io = ((-1).*abs(tend-Ti(index)+c).^(1-p)-(tstart-Ti(index)+c).^(1-p))/(1-p);                
%     elseif (tend-Ti(index)+c)<=0 && (tstart-Ti(index)+c)<=0
%         Io = ((-1).*abs(tend-Ti(index)+c).^(1-p)-(-1).*abs(tstart-Ti(index)+c).^(1-p))/(1-p);        
%     end
end

intLambda = sum(K*Kt*exp(beta*(Mi(index)-Ml)).*Io);
        
end






        




    



















        