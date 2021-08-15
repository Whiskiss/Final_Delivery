%  - - - - - -
%   i a u C p
%  - - - - - -
%
%  Copy a p-vector.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     p             p-vector to be copied
%
%  Returned:
%     c             copy
%
%  This revision:  2008 May 11
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = iauCp(p)

c(1) = p(1);
c(2) = p(2);
c(3) = p(3);

