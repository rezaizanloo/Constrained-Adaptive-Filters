m=5000 ;
n = 4 ;

x=randn(n,m);
W = [3  ; -7 ; 7.5 ; 5];
d=W'*x;

c = [0 ; 0; -1 ; 1.5];
a = 0;

w_lms = zeros(n,m);
w_clms = zeros(n,m);

w_lms = randn(n,1);
w_clms = w_lms;

for j=2:m
    w_clms(:,j) = clms(x(:,j),d(1,j),w_clms(:,j-1),c,a);
    w_lms(:,j) = lms(x(:,j),d(1,j),w_lms(:,j-1),c,a);
end

mse_clms = zeros(n,m); 
mse_lms = zeros(n,m);
close all
for j = 1 : m
mse_clms(:,j)=(w_clms(:,j) - W).^2;
mse_lms(:,j)=(w_lms(:,j) - W).^2;
end
rmse_lms = mean(mse_lms,2);
rmse_clms = mean(mse_clms,2);

sel = 1;  % sel = 1 to 4
t=1:m;
figure,
hold on
plot(t,mse_lms(sel,:),'b');
plot(t,mse_clms(sel,:),'r');
title(['clms by reza izanloo : lms_rmse = ',num2str(rmse_lms(sel)),'clms_rmse = ',num2str(rmse_clms(sel)),]);
xlabel('time');
ylabel('mse');
