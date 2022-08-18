% Joint PDF

function [sampleX,sampleY] = sample_pxy (pxy,x,y)

Fx = cumsum(sum(pxy,1));
        
sampleX = interpola(Fx,x,rand);

indxXUpp = find(x>sampleX,1,'first');
indxXLow = find(x<sampleX,1,'last');

if isempty(indxXLow)
    py_x = pxy(:,1);
elseif isempty(indxXUpp)
    py_x = pxy(:,end);
else
    py_x = interp1([x(indxXLow);x(indxXUpp)],[pxy(:,indxXLow)';pxy(:,indxXUpp)'],sampleX);
end

py_x = py_x/sum(py_x);

Fy_x = cumsum(py_x);    

sampleY = interpola(Fy_x,y,rand);