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
% coco construction
prob = coco_prob();
% prob = coco_set(prob,'coll','NTST',15);
% prob = coco_set(prob,'coll','NCOL',4);

% collocation approximation of ODEs
funcs = {@obv, @obv_dx, @obv_dp, @obv_dxdx, @obv_dxdp, @obv_dpdp};
coll_args = [funcs, {[0; 1], [0 0;0 0], {'l1', 'l2', 'l3'}, [0.0;0.1;0.1]}];
prob = ode_isol2coll(prob, '', coll_args{:});
% boundary conditions
[data, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
bc_funcs = {@obv_bc, @obv_bc_du, @obv_bc_dudu};
prob = coco_add_func(prob, 'bc', bc_funcs{:}, [], 'zero', ...
  'uidx', uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx]));
% objective
data = obj_init_data(data);
prob = coco_add_func(prob, 'obj', @objhan, @objhan_du, data, ...
  'inactive', 'obj', 'uidx', uidx, 'remesh', @obj_remesh);


%% call fmincon
% setup fmincon
options = optimoptions('fmincon','Display','iter');  
options.SpecifyObjectiveGradient  = true;
options.SpecifyConstraintGradient = true;
% options.OptimalityTolerance = 1e-8;
% options.StepTolerance = 1e-8;

u0 = prob.efunc.x0; % Initial point
fprintf('Optimization algorithm: interior-point (default)\n');
x = fmincon(@(u) objfunc(u,prob,'obj'), u0,[],[],[],[],[],[],@(u) nonlincons(u,prob),options);


opt_algorithm = {'sqp', 'sqp-legacy', 'active-set', 'trust-region-reflective'};
for i = 1:numel(opt_algorithm)
    options = optimoptions('fmincon','Display','iter','Algorithm',opt_algorithm{i}); 
    options.SpecifyObjectiveGradient  = true;
    options.SpecifyConstraintGradient = true;
    fprintf('Optimization algorithm: %s\n', opt_algorithm{i});
    x = fmincon(@(u) objfunc(u,prob,'obj'), u0,[],[],[],[],[],[],@(u) nonlincons(u,prob),options);
end



