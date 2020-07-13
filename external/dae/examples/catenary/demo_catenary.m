
% Here we try to solve the problem presented in section 7.3.1 in Recipes.
% We rewrite the x1'=x2 to be x1'=y and y=x1. The Jacobian of
% vector field is computed via numerical differentiation.

t0 = [0; 0.04];
x0 = [1 0; 1 0.04];
y0 = [0; 0.04];
prob = coco_prob();
prob = coco_set(prob,'ddaecoll','Dpoints','Flipped_Gauss_Radau');
prob = coco_set(prob,'ddaecoll','Dnodes','Flipped_Gauss_Radau');
prob = coco_set(prob,'ddaecoll','Apoints','Flipped_Gauss_Radau');
prob = ddaecoll_isol2seg(prob, 'f1', @catenary, t0, x0, y0, []); % Build 'coll' continuation problem
prob = alg_dae_isol2seg(prob, 'f1', @gfunc);

data = coco_get_func_data(prob, 'f1.ddaecoll', 'data'); % Extract toolbox data
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.x1_idx(1); data.T0_idx; data.T_idx], ...
  {'y1s' 'y2s' 'y1e' 'T0' 'T'});
coco(prob, 'coll1', [], 1, {'T' 'y1e'}, [0 1]);

for lab = 1:5
    [sol data] = ddaecoll_read_solution('f1', 'coll1', lab);
    plot(sol.t,sol.x(:,1)); pause(0.5);hold on
end
plot(sol.t,cosh(sol.t),'ro')
figure(10)
plot(sol.t,sinh(sol.t),'ro'); hold on
plot(sol.ta,sol.ya,'b.');

prob = coco_prob();
prob = coco_set(prob, 'ddaecoll', 'NTST', 15);
prob = coco_set(prob,'ddaecoll','Dpoints','Flipped_Gauss_Radau');
prob = coco_set(prob,'ddaecoll','Dnodes','Flipped_Gauss_Radau');
prob = coco_set(prob,'ddaecoll','Apoints','Flipped_Gauss_Radau');
prob = ddaecoll_sol2seg(prob, '', 'coll1', 'f1', 5); % Reconstruct 'coll' continuation problem
prob = alg_dae_sol2seg(prob, '', 'coll1', 'f1', 5);

data = coco_get_func_data(prob, 'ddaecoll', 'data');
% prob = coco_add_func(prob, 'gfunc',@coup, data,'zero','uidx',[data.xbp_idx;data.ybp_idx]);
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.x1_idx(1); data.T0_idx; data.T_idx], ...
  {'y1s' 'y2s' 'y1e' 'T0' 'T'});
coco(prob, 'coll2', [], 1, {'y1e' 'y2s'}, [0 3]);
for lab = 1:11
    figure(2)
    [sol data] = ddaecoll_read_solution('', 'coll2', lab);
    plot(sol.t,sol.x(:,1)); pause(0.5);hold on
end

