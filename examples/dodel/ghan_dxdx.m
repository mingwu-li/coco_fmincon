function dJ = ghan_dxdx(x)

x1 = x(1,:);

dJ = zeros(1,2,2,numel(x1));
dJ(1,1,1,:) = 2;

end