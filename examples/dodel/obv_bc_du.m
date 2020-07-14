function [data, J] = obv_bc_du(prob, data, u) %#ok<INUSD,INUSL>

J = zeros(3,8);
J(1,1) = 1;
J(2,2) = 1;
J(3,4) = 1;

end