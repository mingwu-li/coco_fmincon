%% Stationary points in a two-point boundary-value problem from AUTO [1]
%
% This demo illustrate the use of the wrapper to searching for a local
% extreumum in 1/10*(p1^2+p2^2+p3^2)+\int_0^1 (x1(t)-1)^2 dt to a two-point
% boundary-value problem as follows
%
%     x1' = x2, x2' = -p1*exp(x1+p2*x1^2+p3*x1^4), x1(0) = x1(1) = 0
%
% The optimization problem is first constructed using coco (with coll
% toolbox) and then solved using fmincon.
%
% [1] Doedel, E., Keller, H. B., & Kernevez, J. P. (1991). Numerical
% analysis and control of bifurcation problems (II): Bifurcation in
% infinite dimensions. International Journal of Bifurcation and Chaos,
% 1(04), 745-772.

%% Construct optimization problem
% coco construction
prob = coco_prob();
prob = coco_set(prob,'coll','NTST',20);
prob = coco_set(prob,'coll','NCOL',4);

% collocation approximation of ODEs
funcs = {@obv, @obv_dx, @obv_dp, @obv_dxdx, @obv_dxdp, @obv_dpdp};
coll_args = [funcs, {[0; 1.0], [0 0;0 0], {'p1', 'p2', 'p3'}, [0.0;0.1;0.1]}];
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
  'inactive', 'obj', 'uidx', uidx);


%% Call fmincon
% setup fmincon
options = optimoptions('fmincon','Display','iter');  
options.SpecifyObjectiveGradient  = true;
options.SpecifyConstraintGradient = true;
% options.OptimalityTolerance = 1e-8;
% options.StepTolerance = 1e-8;

u0 = prob.efunc.x0; % Initial point
fprintf('Optimization algorithm: interior-point (default)\n');
x = fmincon(@(u) objfunc(u,prob,'obj'), u0,[],[],[],[],[],[],@(u) nonlincons(u,prob),options);

% different optimization algorithm
opt_algorithm = {'sqp', 'sqp-legacy', 'active-set'};
for i = 1:numel(opt_algorithm)
    options = optimoptions('fmincon','Display','iter','Algorithm',opt_algorithm{i}); 
    options.SpecifyObjectiveGradient  = true;
    options.SpecifyConstraintGradient = true;
    fprintf('Optimization algorithm: %s\n', opt_algorithm{i});
    x = fmincon(@(u) objfunc(u,prob,'obj'), u0,[],[],[],[],[],[],@(u) nonlincons(u,prob),options);
end

%% Post-processing
[~,yy]   = opt_read_sol(x, prob, 'obj');
fprintf('Objetive at located optimum: obj=%d\n', yy);

% optimal state trajectory
coll_sol = opt_read_coll_sol(x, prob, '');
figure(1);
plot(coll_sol.tbp, coll_sol.xbp(:,1), 'r-'); hold on
plot(coll_sol.tbp, coll_sol.xbp(:,2), 'b--');
xlabel('$t$', 'interpreter', 'latex');
legend('$x_1(t)$', '$x_2(t)$', 'interpreter', 'latex');
set(gca, 'Fontsize', 14);

