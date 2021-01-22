function [x] = system_project(x,A,D,d,P,Q)
% Initial estimation error covariance for constrained Kalman filter (system projection)
[u, s, v] = svd(D');
r = length(d); % number of constraints
u2 = u(:, r+1:end);
PND = u2 * u2';
Pc = PND * P * PND;
% Process noise covariance for constrained Kalman filter (system projection).
% Note that this is the *real* process noise covariance.
Qc = PND * Q * PND;

[dQc, lambdaQc, dQcT] = svd(Qc);
 x = A * x + dQc * sqrt(lambdaQc) * randn(size(x));
end