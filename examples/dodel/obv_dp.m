function y = obv_dp(x,p)

x1 = x(1,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
pp = x1+p2.*x1.^2+p3.*x1.^4;
ep = exp(pp);

y = zeros(2,3,numel(x1));
y(2,1,:) = -ep;
y(2,2,:) = -p1.*ep.*x1.^2;
y(2,3,:) = -p1.*ep.*x1.^4;

end