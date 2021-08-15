function [K, x, P] = MeasUpdate(x, z, g, s, G, P, n)

m = length(z);
Inv_W = zeros(m,m);

for i=1:m
    Inv_W(i,i) = s(i)*s(i);    % Inverse weight (measurement covariance)
end

% Kalman gain
K = P*G'*inv(Inv_W+G*P*G');

% State update
x = x + K*(z-g);

% Covariance update
P = (eye(n)-K*G)*P;

