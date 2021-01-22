function [w]=clms(x,d,w,c,a)
% x : input signal
% d : desired signal
% w : unknown parameter
% c & a : constraints : c'*w  = a

u = 1/1000;
y=w'*x;
e=d-y;
w=w+u*x*e;
end