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

%% Construct optimization problem
% Initial trajectory
[t0, x0]  = ode45(@(t,x) linode(t, x, [0.98; 0.3]), [0 2*pi], ...
  [0.276303; 0.960863]);
% coco construction
prob = coco_prob;
prob = coco_set(prob, 'ode', 'autonomous', false);
% collocation approximation of ODEs
coll_args = {@linode, @linode_dx, @linode_dp, @linode_dt,...
    t0, x0, {'k' 'th'}, [0.98; 0.3]};
prob = ode_isol2coll(prob, '', coll_args{:});
% boundary conditions
[data, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
bc_funcs = {@linode_bc, @linode_bc_du};
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
[~,yy]   = opt_read_sol(x, prob, 'vel');
fprintf('Objetive at located optimum: obj=%d\n', yy);

% optimal state trajectory
coll_sol = opt_read_coll_sol(x, prob, '');
figure(1);
plot(coll_sol.tbp, coll_sol.xbp(:,1), 'r-'); hold on
plot(coll_sol.tbp, coll_sol.xbp(:,2), 'b--');
xlabel('$t$', 'interpreter', 'latex');
legend('$x_1(t)$', '$x_2(t)$', 'interpreter', 'latex');
set(gca, 'Fontsize', 14);

figure(2)
plot(coll_sol.xbp(:,1), coll_sol.xbp(:,2), 'ko-', 'MarkerSize', 8);
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'Fontsize', 14);











