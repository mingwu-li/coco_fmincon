function y=obv_dxdp(x,p)

x1 = x(1,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
pp = x1+p2.*x1.^2+p3.*x1.^4;
pd = 1+2*p2.*x1+4*p3.*x1.^3;
ep = exp(pp);

y = zeros(2,2,3,numel(x1));
y(2,1,1,:) = -ep.*pd;
y(2,1,2,:) = -ep.*p1.*x1.*(pd.*x1+2);
y(2,1,3,:) = -ep.*p1.*x1.^3.*(pd.*x1+4);

end