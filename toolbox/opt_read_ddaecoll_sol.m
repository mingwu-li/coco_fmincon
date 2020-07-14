function sol = opt_read_ddaecoll_sol(u, prob, oid)
% This function returns the design variables specific to collocation
% problems with algebraic variables, i.e., x' = f(t,x,y,p) where y can be
% interpreted as control input in optimal control problems. We call x
% differential variables and y algebraic variables.
%
% Here the input u is a vector of all design variables, prob is a problem
% instance of coco, and the string oid should be the object identifyer for
% ddae collocation instance, namely, the id used in ddaecoll_isol2seg(prob, 
% segoid...).

% The output is given by a structure with following fields
%  sol.tbpd - time evaluated at differential base points
%  sol.tbpa - time evaluated at algebraic base points
%  sol.xbp  - states x evaluated at differenaitl base points
%  sol.ybp  - control inputs y evaluted at algebraic base points
%  sol.T0   - initial time
%  sol.T    - period
%  sol.p    - problem parameters

%%
% find the function from prob based on its identifyer
coll_id = coco_get_id(oid, 'ddaecoll');
assert(any(strcmp(coll_id, prob.efunc.identifyers)), 'DDAE collocation id is not found');
Fi = prob.efunc.funcs(strcmp(coll_id, prob.efunc.identifyers));
% extract x from u
x  = u(Fi.x_idx);
% extract data and solution
data = Fi.data;
xbp  = x(data.xbp_idx);
xbp  = reshape(xbp, data.xbp_shp)';
ybp  = x(data.ybp_idx);
ybp  = reshape(ybp, data.y_shp)';
T0   = x(data.T0_idx);
T    = x(data.T_idx);
p    = x(data.p_idx);
tbpd = T0+T*data.tbpd;
tbpa = T0+T*data.tbpa;
% construct sol
sol.tbpd = tbpd;
sol.tbpa = tbpa;
sol.xbp  = xbp;
sol.ybp  = ybp;
sol.T0   = T0;
sol.T    = T;
sol.p    = p;

end