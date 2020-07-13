function J = gfunc_DY(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

% f = x2-y;
J = zeros(1,1,numel(t));
J(1,1,:) = -1;


end