%  - - - - - -
%   i a u I r
%  - - - - - -
%
%  Initialize an r-matrix to the identity matrix.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Returned:
%     r                 r-matrix
%
%  This revision:  2012 April 3
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = iauIr(r)

r(1,1) = 1.0;
r(1,2) = 0.0;
r(1,3) = 0.0;
r(2,1) = 0.0;
r(2,2) = 1.0;
r(2,3) = 0.0;
r(3,1) = 0.0;
r(3,2) = 0.0;
r(3,3) = 1.0;

