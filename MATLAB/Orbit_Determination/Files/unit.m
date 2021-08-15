%--------------------------------------------------------------------------
%
%  unit.m
%
%  this function calculates a unit vector given the original vector. if a
%  zero vector is input, the vector is set to zero.
%
%  input:
%    vec         - vector
%
%  output:
%    outvec      - unit vector
%
%--------------------------------------------------------------------------
function [outvec] = unit ( vec )

small = 0.000001;
magv = norm(vec);

if ( magv > small )
    for i=1:3
        outvec(i)= vec(i)/magv;
    end
else
    for i=1:3
        outvec(i)= 0.0;
    end
end

