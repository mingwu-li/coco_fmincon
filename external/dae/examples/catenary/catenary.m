function f = catenary(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

x1 = x(1,:);
x2 = x(2,:);

f(1,:) = y;
f(2,:) = (1+x2.^2)./x1;

end
