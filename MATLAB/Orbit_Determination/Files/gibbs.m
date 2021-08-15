%--------------------------------------------------------------------------
%
%  gibbs.m
%
%  this function performs the gibbs method of orbit determination. this
%  method determines the velocity at the middle point of the 3 given
%  position vectors.
%
%  inputs:
%    r1          - ijk position vector #1         m
%    r2          - ijk position vector #2         m
%    r3          - ijk position vector #3         m
%
%  outputs:
%    v2          - ijk velocity vector for r2     m/s
%    theta       - angl between vectors           rad
%    error       - flag indicating success        'ok',...
%
%--------------------------------------------------------------------------
function [v2, theta,theta1,copa, error] = gibbs( r1,r2,r3)

global const

small= 0.00000001;
theta= 0.0;
error = '          ok';
theta1= 0.0;

magr1 = norm( r1 );
magr2 = norm( r2 );
magr3 = norm( r3 );
for i= 1:3
    v2(i)= 0.0;
end

p = cross( r2,r3 );
q = cross( r3,r1 );
w = cross( r1,r2 );
pn = unit( p );
r1n = unit( r1 );
copa=  asin( dot( pn,r1n ) );

if ( abs( dot(r1n,pn) ) > 0.017452406 )  
    error= 'not coplanar';
end

d = p + q + w;
magd = norm(d);
n = magr1*p + magr2*q + magr3*w;
magn = norm(n);
nn = unit( n );
dn = unit( d );

% -------------------------------------------------------------
% determine if  the orbit is possible. both d and n must be in
% the same direction, and non-zero.
% -------------------------------------------------------------
if ( ( abs(magd)<small ) || ( abs(magn)<small ) || ...
   ( abs(dot(nn,dn)) < small ) ) % ( dot(nn,dn) < small ) )
    error= '  impossible';
  else
      theta  = angl( r1,r2 );
      theta1 = angl( r2,r3 );

      % ----------- perform gibbs method to find v2 -----------
      r1mr2= magr1-magr2;
      r3mr1= magr3-magr1;
      r2mr3= magr2-magr3;
      s  = r1mr2*r3 + r3mr1*r2 + r2mr3*r1;
      b  = cross( d,r2 );
      l  = sqrt(const.GM_Earth / (magd*magn) );
      tover2= l / magr2;
      v2 = tover2 * b + l * s;
end

