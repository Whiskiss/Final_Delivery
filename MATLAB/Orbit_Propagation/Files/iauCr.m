%  - - - - - -
%   i a u C r
%  - - - - - -
%
%  Copy an r-matrix.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     r            r-matrix to be copied
%
%  Returned:
%   char           copy
%
%  Called:
%     iauCp        copy p-vector
%
%  This revision:  2008 May 11
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = iauCr(r)

c(1,:) = iauCp(r(1,:));
c(2,:) = iauCp(r(2,:));
c(3,:) = iauCp(r(3,:));

