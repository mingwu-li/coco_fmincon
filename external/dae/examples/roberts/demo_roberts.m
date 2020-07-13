
% We solve initial value problem of DAEs presented in 
% https://www.mathworks.com/help/matlab/math/solve-differential-algebraic-equations-daes.html
% y1' = -0.04*y1+10^4*y2*y3
% y2' = 0.04*y1-10^4*y2*y3-3*10^7*y2^2
% 0   = y1+y2+y3-1
% In our framework, we have x=(y1,y2) and y=y3. We introduce p1 and p2 to
% drive them to 10^4 and 10^7 respectively.

% initial solution
y0 = [1.0; 0.0; 0.0];
% tspan = [0 4*logspace(-6,6)];
tspan = 0:0.01:10;
M = [1 0 0; 0 1 0; 0 0 0];
options = odeset('Mass',M,'RelTol',1e-6,'AbsTol',[1e-10 1e-10 1e-10]);
[t,y] = ode15s(@robertsdae,tspan,y0,options);

t0 = t;
x0 = y(:,1:2);
y0 = y(:,3);

prob = coco_prob(); 
prob = coco_set(prob, 'cont', 'h_max', 1e8, 'ItMX', 250);
prob = coco_set(prob, 'ddaecoll', 'NCOL', 4, 'NTST', 50);
prob = coco_set(prob,'ddaecoll','Dpoints','Gauss_Radau');
prob = coco_set(prob,'ddaecoll','Dnodes','Gauss_Radau');
prob = coco_set(prob,'ddaecoll','Apoints','Gauss_Radau');
prob = ddaecoll_isol2seg(prob, '', @roberts, t0, x0, y0, {'p1', 'p2'}, [10 1]); % Build 'impcoll' continuation problem
prob = alg_dae_isol2seg(prob, '', @gfunc);
data = coco_get_func_data(prob, 'ddaecoll', 'data'); % Extract toolbox data
prob = coco_add_pars(prob, 'pars', [data.x0_idx(:); data.T0_idx; data.T_idx], ...
  {'y1s' 'y2s' 'T0' 'T'});
coco(prob, 'impcoll1', [], 1, {'T'}, [4 100]);

[sol, ~] = ddaecoll_read_solution('', 'impcoll1', 5);


y0 = [1.0; 0.0; 0.0];
tspan = 0:0.01:100;
M = [1 0 0; 0 1 0; 0 0 0];
options = odeset('Mass',M,'RelTol',1e-6,'AbsTol',[1e-10 1e-10 1e-10]);
[t,y] = ode15s(@robertsdae,tspan,y0,options);
figure(1);
subplot(1,3,1);
plot(sol.t,sol.x(:,1),'ro'); hold on
plot(t,y(:,1),'r');
subplot(1,3,2)
plot(sol.t,sol.x(:,2),'bo'); hold on
plot(t,y(:,2),'b');
subplot(1,3,3)
plot(sol.t,sol.y(:,1),'ko'); hold on
plot(t,y(:,3),'k');

figure(2)
plot(sol.t,sol.x(:,1)+sol.x(:,2)+sol.y(:,1)-1,'ko'); hold on



% % restart analysis
% prob = coco_prob();
% prob = impcoll_sol2seg(prob, '', 'impcoll1', 4); % Reconstruct 'impcoll' continuation problem
% data = coco_get_func_data(prob, 'impcoll', 'data'); % Extract toolbox data
% prob = coco_add_pars(prob, 'pars', [data.x0idx(:); data.Tidx], ...
%   {'y1s' 'y2s' 'y3s' 'T'});
% coco(prob, 'impcoll2', [], 1, {'p1', 'y3s'}, [1 30]);
% 
% [sol, ~] = impcoll_read_solution('', 'impcoll2', 10);
% figure(1);
% plot(sol.t,sol.x(:,1),'rv'); hold on
% plot(sol.t,sol.x(:,2),'bv');
% plot(sol.t,sol.x(:,3),'kv');
% legend('y1', 'y2', 'y3');
% figure(2)
% plot(sol.t,sol.x(:,1)+sol.x(:,2)+sol.x(:,3)-1,'kv')

