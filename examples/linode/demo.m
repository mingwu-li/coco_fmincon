%% Stationary points in the harmonically excited linear oscillator
%
% This demo illustrate the use of the wrapper to searching for a local
% extreumum in the variable x2(0) along the two-dimensional manifold of
% periodic solutions to the system of differential equations
%
%     x1' = x2, x2' = -x2-k*x1+cos(t+theta)
%
% The optimization problem is first constructed using coco (with coll
% toolbox) and then solved using fmincon.

% The figure shows the one-dimensional solution manifolds obtained in the
% first and third stages of continuation.

%% Construct optimization problem
% Initial trajectory
[t0, x0]  = ode45(@(t,x) linode(t, x, [0.98; 0.3]), [0 2*pi], ...
  [0.276303; 0.960863]);
% coco construction
prob = coco_prob;
prob = coco_set(prob, 'ode', 'autonomous', false);
% collocation approximation of ODEs
coll_args = {@linode, @linode_dx, @linode_dp, @linode_dt, ...
  @linode_dxdx, @linode_dxdp, @linode_dpdp, @linode_dtdx, ...
  @linode_dtdp, @linode_dtdt, t0, x0, {'k' 'th'}, [0.98; 0.3]};
prob = ode_isol2coll(prob, '', coll_args{:});
% boundary conditions
[data, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
bc_funcs = {@linode_bc, @linode_bc_du, @linode_bc_dudu};
prob = coco_add_func(prob, 'po', bc_funcs{:}, [], 'zero', 'uidx', ...
  uidx([maps.x0_idx; maps.x1_idx; maps.T0_idx; maps.T_idx]));
% optimization objective
prob = coco_add_pars(prob, 'vel', uidx(maps.x0_idx(2)), 'v');


%% call fmincon
% setup fmincon
options = optimoptions('fmincon','Display','iter');  
options.SpecifyObjectiveGradient  = true;
options.SpecifyConstraintGradient = true;

u0 = prob.efunc.x0; % initial point
x  = fmincon(@(u) objfunc(u,prob,'vel'), u0,[],[],[],[],[],[],@(u) nonlincons(u,prob),options);

%% results visualization












