% Generate Mi for SEQgen

function Mgen =  generate_M (Ml, Mmax, beta)

Fm = 1-exp(-beta*(Mmax-Ml));

u = rand;

Mgen = -1/beta*log(1-Fm*u)+Ml;

end
