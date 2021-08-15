function [P] = TimeUpdate(P, Phi, Qdt)

if (nargin == 2)
    Qdt = 0.0;
end

P = Phi*P*Phi' + Qdt;