%% Stationary points in the harmonically excited linear oscillator
%
% This demo illustrates the successive continuation approach to searching
% for stationary points that satisfy the necessary conditions for a local
% extremum in the variable x2(0) along the two-dimensional manifold of
% periodic solutions to the system of differential equations 
%
%     x1' = x2, x2' = -x2-k*x1+cos(t+theta)
%
% The augmented continuation problem constructed below has dimensional
% deficit -1-2=-3. The continuation parameters 'k', 'th', 'v', 'd.k',
% 'd.th', and 'd.v' are all initially inactive. Since each stage of
% continuation is along a one-dimensional manifold, four of these
% parameters must be released in each stage. In particular, 'v' must always
% be active. Consequently, 'd.v' should be active during the first two
% stages of continuation, the second of which terminates when this
% parameter equals 1. It should be inactive (and equal to 1) for all
% remaining stages of continuation. Since 'k' is active during all stages
% of continuation, 'd.k' should be inactive during all three stages. Since
% 'th' is inactive during the first two stages of continuation and released
% only in the third stage of continuation, 'd.th' should be active during
% all three stages.
%
% In the first stage of continuation, a local extremum in 'v' is located
% along a one-dimensional solution manifold with trivial Lagrange
% multipliers. This is a branch point, from which emanates a secondary
% one-dimensional submanifold along which the Lagrange multipliers take on
% nontrivial values. As explained above, we terminate continuation along
% this manifold when 'd.v' equals 1. Stationary points within the
% computational domain correspond to points with vanishing 'd.k' and
% 'd.th'.

% The figure shows the one-dimensional solution manifolds obtained in the
% first and third stages of continuation.

%% Initial encoding

prob = coco_prob;
prob = coco_set(prob, 'ode', 'autonomous', false);

%% First run to find local extrema

% zero problems
[t0, x0]  = ode45(@(t,x) linode(t, x, [0.98; 0.3]), [0 2*pi], ...
  [0.276303; 0.960863]);
coll_args = {@linode, @linode_dx, @linode_dp, @linode_dt, ...
  @linode_dxdx, @linode_dxdp, @linode_dpdp, @linode_dtdx, ...
  @linode_dtdp, @linode_dtdt, t0, x0, {'k' 'th'}, [0.98; 0.3]};
prob1 = ode_isol2coll(prob, '', coll_args{:});

[data, uidx] = coco_get_func_data(prob1, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
bc_funcs = {@linode_bc, @linode_bc_du, @linode_bc_dudu};
prob1 = coco_add_func(prob1, 'po', bc_funcs{:}, [], 'zero', 'uidx', ...
  uidx([maps.x0_idx; maps.x1_idx; maps.T0_idx; maps.T_idx]));

prob1 = coco_add_pars(prob1, 'vel', uidx(maps.x0_idx(2)), 'v');


%% fmincon
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp'); 
options = optimoptions('fmincon','Display','iter');  
options.SpecifyObjectiveGradient  = true;
options.SpecifyConstraintGradient = true;

u0 = prob1.efunc.x0;
x = fmincon(@(u) objfunc(u,prob1,'vel'), u0,[],[],[],[],[],[],@(u) nonlincons(u,prob1),options);














