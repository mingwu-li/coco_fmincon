function sol = opt_read_coll_sol(u, prob, oid)
% This function returns the design variables specific to collocation
% problems. 
%
% Here the input u is a vector of all design variables, prob is a problem
% instance of coco, and the string oid should be the object identifyer for
% collocation instance, namely, the id used in ode_isol2coll(prob, segoid..
% .). If bvp toolbox or po toolbox is used, the oid should be in the
% form of ***.bvp.seg1 and ***.po.orb respectively.
%
% The output is given by a structure with following fields
%  sol.tbp - time evaluated at base points
%  sol.xbp - states evaluated at base points
%  sol.T0  - initial time
%  sol.T   - period
%  sol.p   - problem parameters

%%
% find the function from prob based on its identifyer
coll_id = coco_get_id(oid, 'coll');
assert(any(strcmp(coll_id, prob.efunc.identifyers)), 'Collocation id is not found');
Fi = prob.efunc.funcs(strcmp(coll_id, prob.efunc.identifyers));
% extract x from u
x  = u(Fi.x_idx);
% extract data and solution
mesh = Fi.data.coll_seg.mesh;
maps = Fi.data.coll_seg.maps;
tbp = mesh.tbp(maps.tbp_idx);
xbp = x(maps.xbp_idx);
xbp = reshape(xbp, maps.xbp_shp)';
xbp = xbp(maps.tbp_idx,:);
T0  = x(maps.T0_idx);
T   = x(maps.T_idx);
p   = x(maps.p_idx);
% construct sol
sol.tbp = T0+T*tbp;
sol.xbp = xbp;
sol.T0  = T0;
sol.T   = T;
sol.p   = p;

end