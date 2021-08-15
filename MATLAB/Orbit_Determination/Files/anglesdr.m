%---------------------------------------------------------------------------
% 
%  anglesdr.m
%
%  this function solves the problem of orbit determination using three
%  optical sightings.
% 
%  inputs:
%    az1      - azimuth at t1               rad
%    az2      - azimuth at t2               rad
%    az3      - azimuth at t3               rad
%    el1      - elevation at t1             rad
%    el2      - elevation at t2             rad
%    el3      - elevation at t3             rad
%    Mjd1     - Modified julian date of t1
%    Mjd2     - Modified julian date of t2
%    Mjd3     - Modified julian date of t3
%    rsite1   - ijk site1 position vector   m
%    rsite2   - ijk site2 position vector   m
%    rsite3   - ijk site3 position vector   m
%
%  outputs:
%    r        - ijk position vector at t2   m
%    v        - ijk velocity vector at t2   m/s
%
% Last modified:   2015/08/12   M. Mahooti
% 
%---------------------------------------------------------------------------
function [r2,v2] = anglesdr ( az1,az2,az3,el1,el2,el3,Mjd1,Mjd2,Mjd3, ...
                              rsite1,rsite2,rsite3 )

global const eopdata

magr1in = 1.1*const.R_Earth;
magr2in = 1.11*const.R_Earth;
direct  = 'y';

tol    = 1e-8*const.R_Earth;
pctchg = 0.005;

t1 = (Mjd1 - Mjd2)*86400.0;
t3 = (Mjd3 - Mjd2)*86400.0;

los1 = [cos(el1)*sin(az1); cos(el1)*cos(az1); sin(el1)];
los2 = [cos(el2)*sin(az2); cos(el2)*cos(az2); sin(el2)];
los3 = [cos(el3)*sin(az3); cos(el3)*cos(az3); sin(el3)];

[lon1, lat1, h1] = Geodetic(rsite1);
[lon2, lat2, h2] = Geodetic(rsite2);
[lon3, lat3, h3] = Geodetic(rsite3);

M1 = LTC(lon1, lat1);
M2 = LTC(lon2, lat2);
M3 = LTC(lon3, lat3);

% body-fixed system
los1 = M1'*los1;
los2 = M1'*los2;
los3 = M1'*los3;

% mean of date system (J2000)
Mjd_UTC = Mjd1;
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
Mjd_TT = Mjd_UTC + TT_UTC/86400;
Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

P = PrecMatrix(const.MJD_J2000,Mjd_TT);
N = NutMatrix(Mjd_TT);
T = N * P;
E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

los1 = E'*los1;
rsite1 = E'*rsite1;

Mjd_UTC = Mjd2;
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
Mjd_TT = Mjd_UTC + TT_UTC/86400;
Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

P = PrecMatrix(const.MJD_J2000,Mjd_TT);
N = NutMatrix(Mjd_TT);
T = N * P;
E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

los2 = E'*los2;
rsite2 = E'*rsite2;

Mjd_UTC = Mjd3;
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
Mjd_TT = Mjd_UTC + TT_UTC/86400;
Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

P = PrecMatrix(const.MJD_J2000,Mjd_TT);
N = NutMatrix(Mjd_TT);
T = N * P;
E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

los3 = E'*los3;
rsite3 = E'*rsite3;

magr1old  = 99999999.9;
magr2old  = 99999999.9;
magrsite1 = norm(rsite1);
magrsite2 = norm(rsite2);
magrsite3 = norm(rsite3);

cc1 = 2.0*dot(los1,rsite1);
cc2 = 2.0*dot(los2,rsite2);
ktr = 0;

while (abs(magr1in-magr1old) > tol || abs(magr2in-magr2old) > tol)
    ktr = ktr + 1;
    [r2,r3,f1,f2,q1,magr1,magr2,a,deltae32] = doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
                    los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);

    f  = 1.0 - a/magr2*(1.0-cos(deltae32));
    g  = t3 - sqrt(a^3/const.GM_Earth)*(deltae32-sin(deltae32));
    v2 = (r3 - f*r2)/g;
    
    magr1o = magr1in;
    magr1in = (1.0+pctchg)*magr1in;
    deltar1 = pctchg*magr1in;
    [r2,r3,f1delr1,f2delr1,q2,magr1,magr2,a,deltae32] = doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
                           los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);
    pf1pr1 = (f1delr1-f1)/deltar1;
    pf2pr1 = (f2delr1-f2)/deltar1;
 
    magr1in = magr1o;
    deltar1 = pctchg*magr1in;
    magr2o = magr2in;
    magr2in = (1.0+pctchg)*magr2in;
    deltar2 = pctchg*magr2in;
    [r2,r3,f1delr2,f2delr2,q3,magr1,magr2,a,deltae32] = doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
                           los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);
    pf1pr2 = (f1delr2-f1)/deltar2;
    pf2pr2 = (f2delr2-f2)/deltar2;
    
    magr2in = magr2o;
    deltar2 = pctchg*magr2in;
    
    delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
    delta1 = pf2pr2*f1 - pf1pr2*f2;
    delta2 = pf1pr1*f2 - pf2pr1*f1;
    
    deltar1 = -delta1/delta;
    deltar2 = -delta2/delta;
    
    magr1old = magr1in;
    magr2old = magr2in;
    
    magr1in = magr1in + deltar1;
    magr2in = magr2in + deltar2;

end

[r2,r3,f1,f2,q1,magr1,magr2,a,deltae32] = doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
												   los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);

f  = 1.0 - a/magr2*(1.0-cos(deltae32));
g  = t3 - sqrt(a^3/const.GM_Earth)*(deltae32-sin(deltae32));
v2 = (r3 - f*r2)/g;

