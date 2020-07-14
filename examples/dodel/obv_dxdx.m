function y=obv_dxdx(x,p)

x1 = x(1,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
pp = x1+p2.*x1.^2+p3.*x1.^4;
pd = 1+2*p2.*x1+4*p3.*x1.^3;
pt = 2*p2+12*p3.*x1.^2;
ep = exp(pp);

y = zeros(2,2,2,numel(x1));
y(2,1,1,:)=-p1.*ep.*(pt+pd.^2);

end