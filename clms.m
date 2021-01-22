function [w]=clms(x,d,w,c,a)
% x : input signal
% d : desired signal
% w : unknown parameter
% c & a : constraints : c'*w  = a

u = 1/1000;
y=w'*x;
e=d-y;
w=w+2 * u * x * e;
w = w + ((a - c' * w ) / (c' * c)) * c;

end