
% Here we try to solve x''(t)+x'(t)+px(t)=cos(t) with initial conditions.
% In state space, we let x1=x and x2=x', yielding x1'=x2,
% x2'=-x2-p*x1+cos(t). We rewrite the vector field in the following ways
% x1'=x2, x2'=-y-p*x1+cos(t) with algebraic equation y=x1. The Jacobian of
% vector field is provided.


t0 = [0.1; 0.4];
x0 = [1 0.1; 1.2 0.4];
y0 = [0;0.4];

prob = coco_prob();
prob = coco_set(prob,'ddaecoll','Apoints','Gauss');
prob = ddaecoll_isol2seg(prob, '', @linode_het,@linode_het_DFDT,@linode_het_DFDX,@linode_het_DFDY,@linode_het_DFDP, t0, x0, y0,{'p'},1); % Build 'coll' continuation problem
prob = alg_dae_isol2seg(prob, '', @gfunc, @gfunc_DT, @gfunc_DX, @gfunc_DY, @gfunc_DP);
data = coco_get_func_data(prob, 'ddaecoll', 'data'); % Extract toolbox data
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.T0_idx; data.T_idx], ...
  {'y1s' 'y2s' 'T0' 'T'});
coco(prob, 'coll1', [], 1, {'p'}, [0 1.5]);

lab = 4;
[sol data] = ddaecoll_read_solution('', 'coll1', lab);


prob = coco_prob();
prob = coco_set(prob,'ddaecoll', 'Apoints', 'Gauss');
prob = coco_set(prob, 'ddaecoll', 'NTST', 20);
prob = ddaecoll_sol2seg(prob, '', 'coll1', lab);
prob = alg_dae_sol2seg(prob, '', 'coll1', lab);
data = coco_get_func_data(prob, 'ddaecoll', 'data'); % Extract toolbox data
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.T0_idx; data.T_idx], ...
  {'y1s' 'y2s' 'T0' 'T'});
coco(prob, 'coll2', [], 1, {'T'}, [0.1 1.5]);


lab = 4;
[sol data] = ddaecoll_read_solution('', 'coll2', lab);

plot(sol.t,sol.x(:,1))
hold on
plot(sol.t,sol.x(:,2),'r-.')
plot(sol.t,sol.y,'bo')

% symbolic computation 
t0 = 0.1;
p  = 1.5;
syms y(t)
eqn = diff(y,t,2) + diff(y,t) + p*y == cos(t);
Dy = diff(y,t);
cond = [y(t0) == 1,Dy(t0)==0.1];
ySol(t) = dsolve(eqn,cond);
Dysol(t)= diff(ySol,t);

t = 0.1:0.1:1.6;
plot(t,ySol(t),'*');
plot(t,Dysol(t),'ks');