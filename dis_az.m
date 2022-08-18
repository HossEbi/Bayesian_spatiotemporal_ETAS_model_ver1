function [DIS,AZ]=dis_az(LAT1,LON1,LAT2,LON2)
%     -----------------------------------------------------------------
%     This subroutine calculates the distance and azimuth of two points
%     on the earth's surface. An ellipticity correction is made.
%
%     Runstring parameters:
%      LAT1: latitude of first point (negative: south) (f7.3)
%      LON1: longitude of first point (negative: west) (f8.3)
%      LAT2: latitude of second point (negative: south) (f7.3)
%      LON2: longitude of second point (negative: west) (f8.3)
%      DIS:  epicentraldistance [km] (f)
%      ANGI: epicentraldistance [ø] (f)
%      AZ:   azimuth (2nd point relative to 1st point) [ø] (f)
%
%     Reference:
%      H.B. Herrmann (1986)
%      Department of Earth & Atmospheric Sciences
%      221 North Grand Boulevard
%      Saint Louis University
%      St. Louis, Missouri 63103
%
%     Author:
%      P. Smit
%      Imperial College of Science Technology and Medicine
%      Department of Civil and Environmental Engineering
%      Engineering Seismology and Earthquake Engineering Section
%      Imperial College Road
%      London SW7 2BU
%                                                                       
%     History (of FORTRAN program):                                                          
%      18. 3.1991: created by P. Smit                                           
%      31. 7.1998: modified for use at ICSTM by P. Smit     
%     <20. 1.2005: converted to Matlab by J. Douglas
%     -----------------------------------------------------------------

%     ------------------------------------------------------------------
%     initialize some variables

ESQ1= 0.9932315;
RAD= 0.01745329;
EP= LAT1;
EL= LON1;
STNP= LAT2;
STNL= LON2;
SINPR= 0.0;
      
%     ------------------------------------------------------------------
%     eliminate special settings, prepare coordinates

EPIPR= EP* RAD;
if or(EP==90.0,EP==-90.0) 
 ATA= EP* RAD;
elseif or(EP<90.0,EP>-90.0)
 ATA= atan(ESQ1* sin(EPIPR)/ cos(EPIPR));
end

EPIPR= 90.0* RAD- ATA;
EPIPC= EPIPR/ RAD;
if EL<0
 EPILC= 360.0+ EL;
else
 EPILC= EL;
end
EPILR= EPILC* RAD;
STPC= STNP* RAD;
if or(STNP==90.0,STNP==-90.0) 
 ATA= STNP* RAD;
elseif or(STNP<90.0,STNP>-90.0)
 ATA= atan(ESQ1* sin(STPC)/ cos(STPC));
end
STNPR= 90.0* RAD- ATA;
STPC= STNPR/ RAD;
if STNL<0
 STLC= STNL+ 360.0;
else
 STLC= STNL;
end
if STPC<180 
 ANGI= 180.0- EPIPC;
 AZ= 180.0;
 BAZ= 0.0;
else
 STNLR= STLC* RAD;
end

%     ------------------------------------------------------------------
%     convert station coordinates to radians

%     ------------------------------------------------------------------
%     determine polar angle

PANG= abs(STLC- EPILC);
if PANG>180.0 
 PANG= 360.0- PANG;
end
    
PANG= PANG* RAD;
ANGI= cos(EPIPR)* cos(STNPR)+ sin(EPIPR)* sin(STNPR)* cos(PANG);
SNANG=sqrt(abs(1.0- ANGI* ANGI));
ANGI=atan(SNANG/ ANGI);
if ANGI<0
 ANGI= ANGI+ pi;
end
if or(EP==90,EP==-90)
 AZ= STNL* RAD;
end
AZ= (cos(STNPR)- cos(EPIPR)* cos(ANGI))/ (sin(EPIPR)* sin(ANGI));
SNAZ= sqrt(abs(1.0- AZ* AZ));
AZ= atan(SNAZ/ AZ);
if AZ<0
 AZ= AZ+ pi;
end
AZ= AZ* 57.295780;
if and(STNP<90.0,STNP>-90.0) 
 BAZ= (cos(EPIPR)- cos(SINPR)* cos(ANGI))/ (sin(STNPR)* sin(ANGI));
 SNBAZ= sqrt(abs(1.0- BAZ* BAZ));
 BAZ= atan(SNBAZ/ BAZ);
 if BAZ<0
  BAZ= BAZ+ pi;
 end
else
 BAZ= EP* RAD;
end
BAZ= BAZ* 57.295780;
ANGI= ANGI* 57.295780;

%     ------------------------------------------------------------------
%     adjust azimuth and back azimuth

if STLC<EPILC
 if STLC<EPILC-180
  BAZ=360-BAZ;
 elseif STLC==EPILC-180
  if STPC<=180-EPIPC
   AZ=0;
   BAZ=0;
  else
   AZ=180.0;
   BAZ=0.0;
  end
 elseif STPC>EPIPC-180
  AZ= 360.0- AZ;
 end
elseif STLC==EPILC
 if STPC<EPIPC
  AZ=0.0;
  BAZ=180.0;
 elseif STPC==EPIPC
  AZ=0;
  BAZ=0;
 elseif STPC>EPIPC
  AZ=180.0;
  BAZ=0;
 end
elseif STLC>EPILC
 if STLC<EPILC+180
  BAZ= 360.0- BAZ;
 elseif STLC==EPILC+180
  if STPC<=180-EPIPC
   AZ=0;
   BAZ=0;
  else
   AZ=180.0;
   BAZ=0;
  end
 elseif STLC>EPILC+180
  AZ= 360.0- AZ;   
 end
end

%     ------------------------------------------------------------------
%     calculate distance in km

DIS= ANGI* 111.195;
