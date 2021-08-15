%  - - - - - - - -
%   i a u F w 2 m
%  - - - - - - - -
%
%  Form rotation matrix given the Fukushima-Williams angles.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     gamb              F-W angle gamma_bar (radians)
%     phib              F-W angle phi_bar (radians)
%     psi               F-W angle psi (radians)
%     eps               F-W angle epsilon (radians)
%
%  Returned:
%     rotmat            rotation matrix
%
%  Notes:
%  1) Naming the following points:
%
%           e = J2000.0 ecliptic pole,
%           p = GCRS pole,
%           E = ecliptic pole of date,
%     and   P = CIP,
%
%     the four Fukushima-Williams angles are as follows:
%
%        gamb = gamma = epE
%        phib = phi = pE
%        psi = psi = pEP
%        eps = epsilon = EP
%
%  2) The matrix representing the combined effects of frame bias,
%     precession and nutation is:
%
%        NxPxB = R_1(-eps).R_3(-psi).R_1(phib).R_3(gamb)
%
%  3) Three different matrices can be constructed, depending on the
%     supplied angles:
%
%     o  To obtain the nutation x precession x frame bias matrix,
%        generate the four precession angles, generate the nutation
%        components and add them to the psi_bar and epsilon_A angles,
%        and call the present function.
%
%     o  To obtain the precession x frame bias matrix, generate the
%        four precession angles and call the present function.
%
%     o  To obtain the frame bias matrix, generate the four precession
%        angles for date J2000.0 and call the present function.
%
%     The nutation-only and precession-only matrices can if necessary
%     be obtained by combining these three appropriately.
%
%  Called:
%     iauIr        initialize r-matrix to identity
%     iauRz        rotate around Z-axis
%     iauRx        rotate around X-axis
%
%  Reference:
%     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
%
%  This revision:  2009 December 17
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = iauFw2m(gamb, phib, psi, eps)

% Construct the matrix
r = zeros(3);
r = iauIr(r);
r = iauRz(gamb, r);
r = iauRx(phib, r);
r = iauRz(-psi, r);
r = iauRx(-eps, r);

