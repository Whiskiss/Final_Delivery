%  - - - - - - - -
%   i a u P n 0 6
%  - - - - - - - -
%
%  Precession-nutation, IAU 2006 model:  a multi-purpose function,
%  supporting classical (equinox-based) use directly and CIO-based use
%  indirectly.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     date1,date2            TT as a 2-part Julian Date (Note 1)
%     dpsi,deps              nutation (Note 2)
%
%  Returned:
%     epsa                   mean obliquity (Note 3)
%     rb                     frame bias matrix (Note 4)
%     rp                     precession matrix (Note 5)
%     rbp                    bias-precession matrix (Note 6)
%     rn                     nutation matrix (Note 7)
%     rbpn                   GCRS-to-true matrix (Note 8)
%
%  Notes:
%  1)  The TT date date1+date2 is a Julian Date, apportioned in any
%      convenient way between the two arguments.  For example,
%      JD(TT)=2450123.7 could be expressed in any of these ways,
%      among others:
%
%             date1          date2
%
%          2450123.7           0.0       (JD method)
%          2451545.0       -1421.3       (J2000 method)
%          2400000.5       50123.2       (MJD method)
%          2450123.5           0.2       (date & time method)
%
%      The JD method is the most natural and convenient to use in
%      cases where the loss of several decimal digits of resolution
%      is acceptable.  The J2000 method is best matched to the way
%      the argument is handled internally and will deliver the
%      optimum resolution.  The MJD method and the date & time methods
%      are both good compromises between resolution and convenience.
%
%  2)  The caller is responsible for providing the nutation components;
%      they are in longitude and obliquity, in radians and are with
%      respect to the equinox and ecliptic of date.  For high-accuracy
%      applications, free core nutation should be included as well as
%      any other relevant corrections to the position of the CIP.
%
%  3)  The returned mean obliquity is consistent with the IAU 2006
%      precession.
%
%  4)  The matrix rb transforms vectors from GCRS to J2000.0 mean
%      equator and equinox by applying frame bias.
%
%  5)  The matrix rp transforms vectors from J2000.0 mean equator and
%      equinox to mean equator and equinox of date by applying
%      precession.
%
%  6)  The matrix rbp transforms vectors from GCRS to mean equator and
%      equinox of date by applying frame bias then precession.  It is
%      the product rp x rb.
%
%  7)  The matrix rn transforms vectors from mean equator and equinox
%      of date to true equator and equinox of date by applying the
%      nutation (luni-solar + planetary).
%
%  8)  The matrix rbpn transforms vectors from GCRS to true equator and
%      equinox of date.  It is the product rn x rbp, applying frame
%      bias, precession and nutation in that order.
%
%  9)  The X,Y,Z coordinates of the IAU 2000B Celestial Intermediate
%      Pole are elements (3,1-3) of the matrix rbpn.
%
%  10) It is permissible to re-use the same array in the returned
%      arguments.  The arrays are filled in the stated order.
%
%  Called:
%     iauPfw06     bias-precession F-W angles, IAU 2006
%     iauFw2m      F-W angles to r-matrix
%     iauCr        copy r-matrix
%     iauTr        transpose r-matrix
%     iauRxr       product of two r-matrices
%
%  References:
%     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
%     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
%
%  This revision:  2009 December 17
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [epsa, rb, rp, rbp, rn, rbpn] = iauPn06(date1, date2, dpsi, deps)

global const

% Bias-precession Fukushima-Williams angles of J2000.0 = frame bias.
[gamb, phib, psib, eps] = iauPfw06(const.DJM0, const.DJM00);

% B matrix.
r1 = iauFw2m(gamb, phib, psib, eps);
rb = iauCr(r1);

% Bias-precession Fukushima-Williams angles of date.
[gamb, phib, psib, eps] = iauPfw06(date1, date2);

% Bias-precession matrix.
r2 = iauFw2m(gamb, phib, psib, eps);
rbp = iauCr(r2);

% Solve for precession matrix.
rt = iauTr(r1);
rp = iauRxr(r2, rt);

% Equinox-based bias-precession-nutation matrix.
r1 = iauFw2m(gamb, phib, psib + dpsi, eps + deps);
rbpn = iauCr(r1);

% Solve for nutation matrix.
rt = iauTr(r2);
rn = iauRxr(r1, rt);

% Obliquity, mean of date.
epsa = eps;

