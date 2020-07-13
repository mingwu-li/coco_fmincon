function out = robertsdae(t,y)
% p1 = 10^4;
% p2 = 10^7;
p1 = 10;
p2 = 1;
out = [-0.04*y(1) + p1*y(2).*y(3)
   0.04*y(1) - p1*y(2).*y(3) - 3*p2*y(2).^2
   y(1) + y(2) + y(3) - 1 ];