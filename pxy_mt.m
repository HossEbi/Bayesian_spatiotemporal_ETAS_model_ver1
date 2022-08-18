% Joint PDF

function [pxy,sum_intLambda] = pxy_mt (Mi, Ti, ri, m, tstart, tend, Ml, K, theta, dA)

beta  = theta(1);
c     = theta(2);
p     = theta(3);
d     = theta(4);
q     = theta(5);

Krt = (q-1)/pi*d^(2*(q-1))*(p-1)*c^(p-1);

%% Calculate lambda

index = find(Ti < tend);    

lambda = beta*exp(-beta*(m-Ml))*sum((K*exp(beta*(Mi(index)-Ml))./((tend-Ti(index)+c).^p)).*dA.*(Krt./(ri(index).^2+d^2).^q));

%% Calculate intLambda

if p == 1    
    Io = log((tend-Ti(index)+c)./(tstart-Ti(index)+c));     
else    
    Io = ((tend-Ti(index)+c).^(1-p)-(tstart-Ti(index)+c).^(1-p))/(1-p);
end

intLambda = beta*exp(-beta*(m-Ml))*sum((K*exp(beta*(Mi(index)-Ml)).*Io).*dA.*(Krt./(ri(index).^2+d^2).^q));

sum_intLambda = intLambda/(beta*exp(-beta*(m-Ml)));

%% Calculate pxy

pxy = lambda*exp(-intLambda);



        




    



















        