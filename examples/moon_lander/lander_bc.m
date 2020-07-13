function [data, y] = lander_bc(prob, data, u) %#ok<INUSL>

x0 = u(1:2);
x1 = u(3:4);
T0 = u(5);

y = [x0(1)-10; x0(2)+2; x1; T0];

end
