function y = obv_ode(t,x,p1,p2,p3)

x1 = x(1);
x2 = x(2);
y = [x2; -p1*exp(x1+p2*x1.^2+p3*x1.^4)];

end