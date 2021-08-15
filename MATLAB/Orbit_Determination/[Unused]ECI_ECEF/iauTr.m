%  - - - - - -
%   i a u T r
%  - - - - - -
%
%  Transpose an r-matrix.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     r            r-matrix
%
%  Returned:
%     rt           transpose
%
%  Note:
%     It is permissible for r and rt to be the same array.
%
%  Called:
%     iauCr        copy r-matrix
%
%  This revision:  2008 May 22
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rt = iauTr(r)

wm = zeros(3);

for i = 1:3
    for j = 1:3
        wm(i,j) = r(j,i);
    end
end

rt = iauCr(wm);

