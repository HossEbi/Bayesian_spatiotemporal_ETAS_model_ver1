function rxy = calculate_rxy(latitude,longitude,Ggrid)

Grgrid = Ggrid*pi/180; 
rxy = zeros(length(latitude),size(Ggrid,1));

for j=1:length(latitude)      
    Grdata = [latitude(j),longitude(j)]*pi/180;   % Read the aftershock latitude and longitude, gradi a radian
    xv = topgeo(Grgrid(:,1),Grgrid(:,2),Grdata(1),Grdata(2));
    rxy(j,:) = sqrt(xv(:,1).^2+xv(:,2).^2);
%     xj(j,:) = xv(:,1); %mine
%     yj(j,:) = xv(:,2); %mine
end

end