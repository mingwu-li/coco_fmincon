y0 = [0.5; 0.3; 0.2];
% tspan = [0 4*logspace(-6,6)];
tspan = 0:0.01:12;
M = [1 0 0; 0 1 0; 0 0 0];
options = odeset('Mass',M,'RelTol',1e-4,'AbsTol',[1e-6 1e-10 1e-6]);
[t,y] = ode15s(@robertsdae,tspan,y0,options);

figure(1)
plot(t,y(:,1),'r');hold on
plot(t,y(:,2),'b');
plot(t,y(:,3),'k');
legend('y1','y2','y3');
figure(2)
plot(t,y(:,1)+y(:,2)+y(:,3)-1);
