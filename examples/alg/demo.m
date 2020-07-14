% min x^2+y s.t. y=1
clear
%% construct optimization problem
prob = coco_prob;
fcn  = @(p,d,u) deal(d, u(2)-1);      % y-1
obj  = @(p,d,u) deal(d, u(1)^2+u(2)); % x^2+y
% constraint
prob = coco_add_func(prob, 'constraint', fcn, [], 'zero', 'u0', [0;1]);
% objective
prob = coco_add_func(prob, 'obj', obj, [], 'inactive', 'obj', 'uidx', [1;2]);

%% call fmincon
u0 =[1;2]; % initial point
x  = fmincon(@(u) objfunc(u,prob,'obj'), u0,[],[],[],[],[],[],...
    @(u) nonlincons(u,prob));

%% display optimal results
fprintf('Located optimum: x = %d, y = %d\n', x);
% check equality constraint
[~,y1] = opt_read_sol(x,prob,'constraint');
fprintf('Evaluation of constraint at located optimum: h=%d\n', y1);
% optimal objective
[~,y2] = opt_read_sol(x,prob,'obj');
fprintf('Objetive at located optimum: obj=%d\n', y2);
