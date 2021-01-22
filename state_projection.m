% constrained kalman filter
% atate Projectionj method (W=INV(P))
% By Reaz izanloo 
close all
clear all
clc

tf = 300; % final time (seconds)
T = 3; % time step (seconds)
theta = pi / 3; % heading angle (measured CCW from east)
tantheta = tan(theta);
sintheta = sin(theta);
costheta = cos(theta);

% System matrix.
A = [1 0 T 0
   0 1 0 T
   0 0 1 0
   0 0 0 1];

% State constraint matrix.
D = [1 -tantheta 0 0; 
        0 0 1 -tantheta];
d = [0;0] ;

% new matrix B
 B =[1 0 0 0
    0 1 0 0
    0 0 0 1];

Q = diag([4, 4, 1, 1]); % Process noise covariance (m, m, m/sec, m/sec)
Qsqrt = sqrt(Q);

R = diag([3, 3,1]); % Measurement noise covariance (m, m)
Rsqrt = sqrt(R);

x = [0; 0; tantheta; 1] * 100;
xhat=x;
xhat_ckf= x;
xhat_kf=x;
% Initial estimation error covariance
p_ckf = diag([R(1,1), R(2,2), Q(1,1), Q(2,2)]);
p_kf=p_ckf;

xhat_kf_sequence=[];
xhat_ckf_sequence=[];
xarray=x;

% Begin the simulation.
for t = T : T : tf
 
   % Get the measurement.
   z=B*x;
   MeasErr = Rsqrt*randn(size(z));
   z = z + MeasErr;
     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run kf
   xhat_kf=A*xhat_kf;
   p_kf=A*p_kf*A'+Q;
   
   K=p_kf*B'*inv(B*p_kf*B'+R);
   xhat_kf=xhat_kf+K*(z-B*xhat_kf);
   p_kf=(eye(4)-K*B)*p_kf;
   
   xhat_kf_sequence=[xhat_kf_sequence xhat_kf];
   
  %%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%% run ckf 
     xhat_ckf=A*xhat_ckf;
  p_ckf=A*p_ckf*A'+Q;
  
   K=p_ckf*B'*inv(B*p_ckf*B'+R);
   xhat_ckf=xhat_ckf+K*(z-B*xhat_ckf);
     p_ckf=(eye(4)-K*B)*p_ckf;
      
   xhat_ckf=xhat_ckf-p_ckf*D'*inv(D*p_ckf*D')*D*xhat_ckf;
%    p_ckf = (eye(4) -p_ckf*D'*inv(D*p_ckf*D')*D) * p_ckf * (eye(4) -p_ckf*D'*inv(D*p_ckf*D')*D)';

 xhat_ckf_sequence=[xhat_ckf_sequence xhat_ckf];
    
  % Simulate the system.
   x = A*x+ Qsqrt*randn(size(x));
   % Constrain the vehicle (i.e., the true state) to the straight road.
   if abs(x(1) - tantheta * x(2)) > 2
      x(2) = (x(2) + x(1) * tantheta) / (1 + tantheta^2);
      x(1) = x(2) * tantheta;
   end
   if abs(x(3) - tantheta * x(4)) > 0.2
      x(4) = (x(4) + x(3) * tantheta) / (1 + tantheta^2);
      x(3) = x(4) * tantheta;
   end
   xarray = [xarray x];
   
end
xarray(:,floor(tf/T)+1)=[];

% Plot data.

c=1:4;
t = 0 : T : tf-T;
se_kf=(xarray- xhat_kf_sequence).^2;
RMSE=sqrt(mean(se_kf,2));
se_ckf=(xarray- xhat_ckf_sequence).^2;
RMSE2=sqrt(mean(se_ckf,2));

sel = 3;
figure,
hold on
plot(t,sqrt(se_kf(sel,:)) ,'K:');
plot(t,sqrt(se_ckf(sel,:)) ,'K-');
legend('kf','ckf');
title([' RMSE of estimation :    Unconstrained,RMSE = ',num2str(RMSE(sel)),'     constrained, RMSE = ', num2str(RMSE2(sel)),'']);
xlabel('seconds');
ylabel('RMSE');

