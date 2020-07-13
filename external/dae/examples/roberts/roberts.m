function f = roberts(t,x,y,p)
%CATENARx   'coll'-compatible encoding of catenarx vector field

x1 = x(1,:);
x2 = x(2,:);
y1 = y(1,:);
p1 = p(1,:);
p2 = p(2,:);

f(1,:) = -0.04*x1 + p1.*x2.*y1;
f(2,:) = 0.04*x1 - p1.*x2.*y1 - 3*p2.*x2.^2;

end