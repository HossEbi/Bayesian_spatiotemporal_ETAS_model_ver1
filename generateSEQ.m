% Function for generating SEQgen

function [Mseq,timeseq,Lonseq,Latseq,rseq,sum_intLambdaseq] = generateSEQ (Mi, Ti, ri, tstart, tend, Ml, sampleK, sampletheta, Xgrid, Ygrid, Ggrid, Mmax, dA)
%Mi=M(indexSeq); Ti=time_MS(indexSeq);ri=rxy(indexSeq,:);Ml=Mc;sampleK=K(j);sampletheta=samples(:,j);Xgrid=Xcgrid; Ygrid=Ycgrid;

%% Generating 1st Magnitude, time, and Location

Mseq    = [];
timeseq = [];
Lonseq  = [];
Latseq  = []; 
rseq    = []; 
sum_intLambdaseq = [];

[Mgen,tgen,numRepeat] = generate_m_time (Mi, Ti, ri, tstart, tend, Ml, sampleK, sampletheta, Mmax);

if tgen > tend
    
    display(['                  tgen = ',num2str(tgen),' > T_end = ',num2str(tend),', NO sequence will be generated!!!'])
    
else
    
    [rgen,Longen,Latgen,sum_intLambda] = generate_R(Mi, Ti, ri, Mgen, tstart, tgen, Ml, sampleK, sampletheta, Xgrid, Ygrid, Ggrid, dA);

    count = 0;
    while tgen <= tend
        
        display(['                  tgen = ',num2str(tgen),' / ',num2str(tend),' , Mgen = ',num2str(Mgen),' - Thinning iteration = ',num2str(numRepeat)])
   
        count = count + 1;
        Mseq(count,:)    = Mgen;
        timeseq(count,:) = tgen;
        Lonseq(count,:)  = Longen;
        Latseq(count,:)  = Latgen;
        rseq(count,:)    = rgen;
        sum_intLambdaseq(count,:) = sum_intLambda;
    
        [Mgen,tgen,numRepeat] = generate_m_time ([Mi;Mseq],[Ti;timeseq],[ri;rseq], tgen, tend, Ml, sampleK, sampletheta, Mmax);
        [rgen,Longen,Latgen,sum_intLambda] = generate_R([Mi;Mseq],[Ti;timeseq],[ri;rseq], Mgen, timeseq(end), tgen, Ml, sampleK, sampletheta, Xgrid, Ygrid, Ggrid, dA);
        
    end

end
 
end