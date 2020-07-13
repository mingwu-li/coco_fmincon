function y = lander0(t,x,p)

y = [x(2,:); -1.5+t.*p(1,:)];

end