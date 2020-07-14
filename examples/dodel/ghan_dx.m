function J = ghan_dx(x)

x1 = x(1,:);

J  = zeros(1,2,numel(x1));
J(1,1,:) = 2*(x1-1);

end
