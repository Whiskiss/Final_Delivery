function [r2,r3,f1,f2,q1,magr1,magr2,a,deltae32] = doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
                      los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct)

global const

rho1 = (-cc1 + sqrt(cc1^2-4*(magrsite1^2-magr1in^2))) / 2.0;
rho2 = (-cc2 + sqrt(cc2^2-4*(magrsite2^2-magr2in^2))) / 2.0;

r1 = rho1.*los1 + rsite1;
r2 = rho2.*los2 + rsite2;

magr1 = norm(r1);
magr2 = norm(r2);
 
if direct == 'y'
    w = cross(r1,r2)/(magr1*magr2);
else
    w = -cross(r1,r2)/(magr1*magr2);
end

rho3 =  -dot(rsite3,w)/dot(los3,w);
r3 = rho3.*los3 + rsite3;
magr3 = norm(r3);

cosdv21 = dot(r2,r1)/(magr2*magr1);
sindv21 = norm(cross(r2,r1))/(magr2*magr1);
dv21 = atan2(sindv21,cosdv21);

cosdv31 = dot(r3,r1)/(magr3*magr1);
sindv31 = sqrt(1.0 - cosdv31^2);
dv31 = atan2(sindv31,cosdv31);

cosdv32 = dot(r3,r2)/(magr3*magr2);
sindv32 = norm(cross(r3,r2))/(magr3*magr2);
dv32 = atan2(sindv32,cosdv32);    

if dv31 > pi
    c1 = (magr2*sindv32)/(magr1*sindv31);
    c3 = (magr2*sindv21)/(magr3*sindv31);
    p = (c1*magr1+c3*magr3-magr2)/(c1+c3-1);
else
    c1 = (magr1*sindv31)/(magr2*sindv32);
    c3 = (magr1*sindv21)/(magr3*sindv32);
    p = (c3*magr3-c1*magr2+magr1)/(-c1+c3+1);
end

ecosv1 = p/magr1-1;
ecosv2 = p/magr2-1;
ecosv3 = p/magr3-1;

if dv21 ~= pi
    esinv2 = (-cosdv21*ecosv2+ecosv1)/sindv21;
else
    esinv2 = (cosdv32*ecosv2-ecosv3)/sindv31;
end

e = sqrt(ecosv2^2+esinv2^2);
a = p/(1-e^2);

if e*e < 0.99
    n = sqrt(const.GM_Earth/a^3);
    
    s = magr2/p*sqrt(1-e^2)*esinv2;
    c = magr2/p*(e^2+ecosv2);
    
    sinde32 = magr3/sqrt(a*p)*sindv32-magr3/p*(1-cosdv32)*s;
    cosde32 = 1-magr2*magr3/(a*p)*(1-cosdv32);
    deltae32 = atan2(sinde32,cosde32);
    
    sinde21 = magr1/sqrt(a*p)*sindv21+magr1/p*(1-cosdv21)*s;
    cosde21 = 1-magr2*magr1/(a*p)*(1-cosdv21);
    deltae21 = atan2(sinde21,cosde21);
    
    deltam32 = deltae32+2*s*(sin(deltae32/2))^2-c*sin(deltae32);
    deltam12 = -deltae21+2*s*(sin(deltae21/2))^2+c*sin(deltae21);
else
    n = sqrt(GM_Earth/-a^3);
    
    s = magr2/p*sqrt(e^2-1)*esinv2;
    c = magr2/p*(e^2+ecosv2);
    
    sindh32 = magr3/sqrt(-a*p)*sindv32-magr3/p*(1-cosdv32)*s;
    sindh21 = magr1/sqrt(-a*p)*sindv21+magr1/p*(1-cosdv21)*s;
    
    deltah32 = log( sindh32 + sqrt(sindh32^2 +1) );
    deltah21 = log( sindh21 + sqrt(sindh21^2 +1) );
    
    deltam32 = -deltah32+2*s*(sinh(deltah32/2))^2+c*sinh(deltah32);
    deltam12 = deltah21+2*s*(sinh(deltah21/2))^2-c*sinh(deltah21);
    
    deltae32=deltah32;
end

f1 = t1-deltam12/n;
f2 = t3-deltam32/n;

q1 = sqrt(f1^2+f2^2);

