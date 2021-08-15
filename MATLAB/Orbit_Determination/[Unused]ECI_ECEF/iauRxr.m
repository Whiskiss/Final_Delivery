%  - - - - - - -
%   i a u R x r
%  - - - - - - -
%
%  Multiply two r-matrices.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     a            first r-matrix
%     b            second r-matrix
%
%  Returned:
%     atb          a * b
%
%  Note:
%     It is permissible to re-use the same array for any of the
%     arguments.
%
%  Called:
%     iauCr        copy r-matrix
%
%  This revision:  2008 November 18
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function atb = iauRxr(a, b)

wm = zeros(3);

for i = 1:3
    for j = 1:3
        w = 0;
        for k = 1:3
            w = w + a(i,k) * b(k,j);
        end
        wm(i,j) = w;
    end
end

atb = iauCr(wm);

