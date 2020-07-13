

% moon lander optimal control problem
% minimize \int_0^tf u(t)dt
% s.t.     \dot{x}_1 = x_2, \dot{x}_2 = -1.5 + u(t), 0<=u(t)<=3
%          x_1(0) = 10, x_2(0)=2, x_1(tf)=0, x_2(tf)=0
% Here tf is free

%% construct initial solution with free tf and u = k*t
t0 = 0;
x0 = [10 -2];
prob = coco_prob();
prob = coco_set(prob, 'cont', 'h_max', 100);
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = ode_isol2coll(prob, '', @lander0, t0, x0, 'up', 0);
data = coco_get_func_data(prob, 'coll', 'data'); % Extract function data
maps = data.coll_seg.maps;
prob = coco_add_pars(prob, 'pars', ...
  [maps.x0_idx; maps.x1_idx; maps.T0_idx; maps.T_idx], ...
  {'x1s' 'x2s' 'x1e' 'x2e' 'T0' 'T'});
cont_args = {1, {'x1e' 'T' 'x2e'}, {[0 20],[0, 100]}};
fprintf('\n Run=''%s'': Continue trajectory segments until x1e=0.\n', ...
  'coll1');
bd = coco(prob, 'coll1', [], cont_args{:});

%
labs = coco_bd_labs(bd, 'EP');
prob = coco_prob();
prob = coco_set(prob, 'cont', 'h_max', 100);
prob = ode_coll2coll(prob, '', 'coll1', max(labs));
data = coco_get_func_data(prob, 'coll', 'data');
maps = data.coll_seg.maps;
prob = coco_add_pars(prob, 'pars', ...
  [maps.x0_idx; maps.x1_idx; maps.T0_idx; maps.T_idx], ...
  {'x1s' 'x2s' 'x1e' 'x2e' 'T0' 'T'});
cont_args = {1, {'x2e' 'up' 'T'}, {[-6 0],[-100, 10]}};
fprintf(...
  '\n Run=''%s'': Continue segments from point %d in run ''%s'' until x2e=0.\n', ...
  'coll2', max(labs), 'coll1');
bd = coco(prob, 'coll2', [], cont_args{:});


%% construct optimal control problem
labs = coco_bd_labs(bd, 'EP');
sol  = coll_read_solution('', 'coll2', max(labs));
t0   = sol.tbp;
x0   = sol.xbp;
y0   = sol.p*sol.tbp;
prob = coco_prob();
prob = coco_set(prob,'ddaecoll','Apoints','Gauss');
prob = ddaecoll_isol2seg(prob, '', @lander, t0, x0, y0, []); 
prob = alg_dae_isol2seg(prob, 'g1', '', @g1func, 'inequality'); % u>=0
prob = alg_dae_isol2seg(prob, 'g2', '', @g2func, 'inequality'); % u<=3

[data, uidx] = coco_get_func_data(prob, 'ddaecoll', 'data', 'uidx');
bc_funcs = {@lander_bc, @lander_bc_du};
prob = coco_add_func(prob, 'bc', bc_funcs{:}, [], 'zero', 'uidx', ...
  uidx([data.x0_idx; data.x1_idx; data.T0_idx]));

prob = coco_add_func(prob, 'obj', @int_u, data, 'inactive', 'obj', 'uidx', uidx([data.ybp_idx; data.T_idx]));


%% fmincon
options = optimoptions('fmincon','Display','iter');  
options.SpecifyObjectiveGradient  = true;
options.SpecifyConstraintGradient = true;
% options.OptimalityTolerance = 1e-8;
% options.StepTolerance = 1e-8;

u0 = prob.efunc.x0;
fprintf('Optimization algorithm: interior-point (default)\n');
x = fmincon(@(u) objfunc(u,prob,'obj'), u0,[],[],[],[],[],[],@(u) nonlincons(u,prob),options);

xx = x(data.xbp_idx);
xx = reshape(xx, data.xbp_shp);
yy = x(data.ybp_idx);
T0 = x(data.T0_idx);
T  = x(data.T_idx);
tbpd = T0+T*data.tbpd;
tbpa = T0+T*data.tbpa;

figure(1)
plot(tbpd, xx(1,:)); hold on
plot(tbpd, xx(2,:));

figure(2)
plot(tbpa, yy);

% opt_algorithm = {'sqp', 'sqp-legacy', 'active-set', 'trust-region-reflective'};
% for i = 1:numel(opt_algorithm)
%     options = optimoptions('fmincon','Display','iter','Algorithm',opt_algorithm{i}); 
%     options.SpecifyObjectiveGradient  = true;
%     options.SpecifyConstraintGradient = true;
%     fprintf('Optimization algorithm: %s\n', opt_algorithm{i});
%     x = fmincon(@(u) objfunc(u,prob1,'obj'), u0,[],[],[],[],[],[],@(u) nonlincons(u,prob1),options);
% end










