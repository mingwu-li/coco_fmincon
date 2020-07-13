% what we need
% min x^2+y s.t. y=1

prob = coco_prob;
fcn  = @(p,d,u) deal(d, u(2)-1);
obj  = @(p,d,u) deal(d, u(1)^2+u(2));

prob = coco_add_func(prob, 'nc', fcn, [], 'zero', 'u0', [0;1]);
prob = coco_add_func(prob, 'obj', obj, [], 'inactive', 'obj', 'uidx', [1;2]);


% [x,y]=[1,2];

x = fmincon(@(u) objfunc(u,prob,'obj'), [1;2],[],[],[],[],[],[],@(u) nonlincons(u,prob));







