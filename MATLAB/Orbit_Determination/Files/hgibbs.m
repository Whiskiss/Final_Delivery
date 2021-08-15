%--------------------------------------------------------------------------
%
%  hgibbs.m
%
%  this function implements the herrick-gibbs approximation for orbit
%  determination, and finds the middle velocity vector for the 3 given
%  position vectors.
%
%  inputs:
%    r1          - ijk position vector #1         m
%    r2          - ijk position vector #2         m
%    r3          - ijk position vector #3         m
%    Mjd1        - julian date of 1st sighting    days from 4713 bc
%    Mjd2        - julian date of 2nd sighting    days from 4713 bc
%    Mjd3        - julian date of 3rd sighting    days from 4713 bc
%
%  outputs:
%    v2          - ijk velocity vector for r2     m/s
%    theta       - angl between vectors           rad
%    error       - flag indicating success        'ok',...
%
%--------------------------------------------------------------------------
function [v2, theta,theta1,copa, error] = hgibbs (r1,r2,r3,Mjd1,Mjd2,Mjd3)

SAT_Const

error =  '          ok';
theta = 0.0;
theta1= 0.0;
magr1 = norm( r1 );
magr2 = norm( r2 );
magr3 = norm( r3 );

for i= 1 : 3
    v2(i)= 0.0;
end

tolangle= 0.01745329251994;
dt21= (Mjd2-Mjd1)*86400.0;
dt31= (Mjd3-Mjd1)*86400.0;
dt32= (Mjd3-Mjd2)*86400.0;

p = cross( r2,r3 );
pn = unit( p );
r1n = unit( r1 );
copa=  asin( dot( pn,r1n ) );

if ( abs( dot(r1n,pn) ) > 0.017452406 )
    error= 'not coplanar';
end

theta  = angl( r1,r2 );
theta1 = angl( r2,r3 );

if ( (theta > tolangle) | (theta1 > tolangle) )  
    error= '   angl > 1ø';
end

term1= -dt32*( 1.0/(dt21*dt31) + GM_Earth/(12.0*magr1*magr1*magr1) );
term2= (dt32-dt21)*( 1.0/(dt21*dt32) + GM_Earth/(12.0*magr2*magr2*magr2) );
term3=  dt21*( 1.0/(dt32*dt31) + GM_Earth/(12.0*magr3*magr3*magr3) );

v2 =  term1*r1 + term2* r2 + term3* r3;

