function J = gfunc_DX(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

% f = x2-y;
J = zeros(1,2,numel(t));
J(1,2,:) = 1;


end