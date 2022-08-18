function [rgen,Longen,Latgen,sum_intLambda] = generate_R(Mi, Ti, ri, Mgen, tstart, tgen, Ml, K, theta, Xgrid, Ygrid, Ggrid, dA)
% ri=rxy
pxy = zeros(length(Ggrid),1);
sum_intLambda = zeros(1,length(Ggrid));
for j=1:length(Ggrid)
    [pxy(j),sum_intLambda(j)] = pxy_mt (Mi, Ti, ri(:,j), Mgen, tstart, tgen, Ml, K, theta, dA);
end
pxy = (reshape(pxy,length(Xgrid),length(Ygrid)))'/sum(pxy);

[Longen,Latgen] = sample_pxy(pxy,Xgrid,Ygrid);

rgen = calculate_rxy(Latgen,Longen,Ggrid);

end

