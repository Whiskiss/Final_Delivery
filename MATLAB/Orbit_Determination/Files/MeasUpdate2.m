function [K, x, P] = MeasUpdate2(x, z, g, s, G, P, n, delta, K_prev, lim)

m = length(z);
Inv_W = zeros(m,m);

for i=1:m
    Inv_W(i,i) = s(i)*s(i);    % Inverse weight (measurement covariance)
end

% Kalman gain
if delta>lim
    K = K_prev;
else
    K = P*G'*inv(Inv_W+G*P*G');
end

% State update
x = x + K*(z-g);

% Covariance update
P = (eye(n)-K*G)*P;

