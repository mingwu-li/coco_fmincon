function demo_adjoint()
clear
prob = coco_prob();
prob = coco_set(prob,'cont', 'NAdapt', 5);
p1 = [-0.1, 3.5];
p2 = [-0.1, 0.8];
p3 = [-0.1, 0.6];
objint = [0, 1.5];


%% first run to find initial fold
% zero problems
funcs = {@obv, @obv_dx, @obv_dp, @obv_dxdx, @obv_dxdp, @obv_dpdp};
coll_args = [funcs, {[0; 1], [0 0;0 0], {'l1', 'l2', 'l3'}, [0.0;0.1;0.1]}];
prob1 = ode_isol2coll(prob, '', coll_args{:});

[data, uidx] = coco_get_func_data(prob1, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
bc_funcs = {@obv_bc, @obv_bc_du, @obv_bc_dudu};
prob1 = coco_add_func(prob1, 'bc', bc_funcs{:}, [], 'zero', ...
  'uidx', uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx]));

data  = obj_init_data(data);
prob1 = coco_add_func(prob1, 'obj', @objhan, @objhan_du, data, ...
  'inactive', 'obj', 'uidx', uidx, 'remesh', @obj_remesh);

% adjoints
prob1 = adjt_isol2coll(prob1, '');

[data, axidx] = coco_get_adjt_data(prob1, 'coll', 'data', 'axidx');
opt = data.coll_opt;
prob1 = coco_add_adjt(prob1, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx]));

data  = adj_obj_init_data(data);
prob1 = coco_add_adjt(prob1, 'obj', @adj_objhan, @adj_objhan_du, data, ...
  'd.obj', 'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]), ...
  'remesh', @adj_obj_remesh);

cont_pars =  {'obj', 'l1', 'd.obj', 'd.coll.T0', 'd.l3', 'd.l2', 'd.l1'};
coco(prob1, 'obv1', [], 1, cont_pars, {objint,p1});

%% branch switch from fold to grow nontrivial adjoint
bd1   = coco_bd_read('obv1');
BPlab = coco_bd_labs(bd1, 'BP');
lab   = BPlab(3);

% zero problems
prob2 = ode_BP2coll(prob, '', 'obv1', lab);

[data, uidx] = coco_get_func_data(prob2, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob2 = coco_add_func(prob2, 'bc', bc_funcs{:}, [], 'zero', ...
  'uidx', uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx]));

data  = obj_init_data(data);
prob2 = coco_add_func(prob2, 'obj', @objhan, @objhan_du, data, ...
  'inactive', 'obj', 'uidx', uidx, 'remesh', @obj_remesh);

% branch switch data
chart = coco_read_solution('obv1', lab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');

% adjoints
prob2 = adjt_BP2coll(prob2, '', 'obv1', lab);

[chart, aidx] = coco_read_adjoint('bc', 'obv1', lab, 'chart', 'lidx');
[data, axidx] = coco_get_adjt_data(prob2, 'coll', 'data', 'axidx');
opt = data.coll_opt;
prob2   = coco_add_adjt(prob2, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx]), ...
  'l0', chart.x, 'tl0', cdata.v(aidx));

[chart, aidx] = coco_read_adjoint('obj', 'obv1', lab, 'chart', 'lidx');
data  = adj_obj_init_data(data);
prob2 = coco_add_adjt(prob2, 'obj', @adj_objhan, @adj_objhan_du, data, ...
  'd.obj', 'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]),...
  'l0', chart.x, 'tl0', cdata.v(aidx), 'remesh', @adj_obj_remesh);

% computational domain
dobj_int = [chart.x(1) 1.1];
prob2 = coco_add_event(prob2, 'opt', 'BP', 'd.obj', '>', 1);

cont_pars = {'d.obj', 'obj', 'l1', 'd.coll.T0', 'd.l3', 'd.l2', 'd.l1'};
% cont_pars = {'d.obj', 'obj', 'l2', 'd.coll.T0', 'd.l1', 'd.l3', 'd.l2'};
% cont_pars = {'d.obj', 'obj', 'l3', 'd.coll.T0', 'd.l1', 'd.l2', 'd.l3'};
coco(prob2, 'obv2', [], cont_pars, dobj_int);

%% continue to let d.l2=0
bd2 = coco_bd_read('obv2');
lab = coco_bd_labs(bd2, 'opt');

% zero problems
prob3 = ode_coll2coll(prob, '', 'obv2', lab);

[data, uidx] = coco_get_func_data(prob3, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob3 = coco_add_func(prob3, 'bc', bc_funcs{:}, [], 'zero', 'uidx', ...
  uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx]));

data  = obj_init_data(data);
prob3 = coco_add_func(prob3, 'obj', @objhan, @objhan_du, data, ...
  'inactive', 'obj', 'uidx', uidx, 'remesh', @obj_remesh);

% adjoints
prob3 = adjt_coll2coll(prob3, '', 'obv2', lab);

chart = coco_read_adjoint('bc', 'obv2', lab, 'chart');
[data, axidx] = coco_get_adjt_data(prob3, 'coll', 'data', 'axidx');
opt = data.coll_opt;
prob3   = coco_add_adjt(prob3, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx]), ...
  'l0', chart.x);

chart = coco_read_adjoint('obj', 'obv2', lab, 'chart');
data  = adj_obj_init_data(data);
prob3 = coco_add_adjt(prob3, 'obj', @adj_objhan, @adj_objhan_du, data,...
  'd.obj', 'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]),...
  'l0', chart.x, 'remesh', @adj_obj_remesh);

prob3 = coco_add_event(prob3, 'opt', 'd.l2', '=', 0);
cont_pars = {'l2', 'l1', 'obj', 'd.l2', 'd.l3', 'd.coll.T0', 'd.l1'};
coco(prob3, 'obv3', [], cont_pars, {p2, p1, objint});

%% continue to let d.l3=0
bd3 = coco_bd_read('obv3');
lab = coco_bd_labs(bd3, 'opt');
lab = lab(1);
% zero problems
prob4 = ode_coll2coll(prob, '', 'obv3', lab);

[data, uidx] = coco_get_func_data(prob4, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob4 = coco_add_func(prob4, 'bc', bc_funcs{:}, [], 'zero', 'uidx', ...
  uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx]));

data  = obj_init_data(data);
prob4 = coco_add_func(prob4, 'obj', @objhan, @objhan_du, data, ...
  'inactive', 'obj', 'uidx', uidx, 'remesh', @obj_remesh);

% adjoints
prob4 = adjt_BP2coll(prob4, '', 'obv3', lab);

chart = coco_read_adjoint('bc', 'obv3', lab, 'chart');
[data, axidx] = coco_get_adjt_data(prob4, 'coll', 'data', 'axidx');
opt = data.coll_opt;
prob4   = coco_add_adjt(prob4, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx]), ...
  'l0', chart.x);

chart = coco_read_adjoint('obj', 'obv3', lab, 'chart');
data  = adj_obj_init_data(data);
prob4 = coco_add_adjt(prob4, 'obj', @adj_objhan, @adj_objhan_du, data,...
  'd.obj', 'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]),...
  'l0', chart.x, 'remesh', @adj_obj_remesh);

prob4 = coco_add_event(prob4, 'opt', 'd.l3', '=', 0);
cont_pars = {'l3', 'l2', 'l1', 'obj', 'd.l3', 'd.coll.T0', 'd.l1'};
coco(prob4, 'obv4', [], cont_pars, {p3,p2,p1,objint});


%% plot results
figure(1)
bd  = coco_bd_read('obv4');
lab = coco_bd_labs(bd, 'opt');
sol = coll_read_solution('', 'obv4', lab);
plot(sol.tbp, sol.xbp(:,1), 'ro'); hold on
plot(sol.tbp, sol.xbp(:,2), 'bv');
grid on 
box on

end

function y = obv(x,p)

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);

y(1,:) = x2;
y(2,:) = -p1.*exp(x1+p2.*x1.^2+p3.*x1.^4);

end

function y = obv_dx(x,p)

x1 = x(1,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
pp = x1+p2.*x1.^2+p3.*x1.^4;
pd = 1+2*p2.*x1+4*p3.*x1.^3;

y = zeros(2,2,numel(x1));
y(1,2,:) = 1;
y(2,1,:) = -p1.*exp(pp).*pd;

end

function y = obv_dp(x,p)

x1 = x(1,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
pp = x1+p2.*x1.^2+p3.*x1.^4;
ep = exp(pp);

y = zeros(2,3,numel(x1));
y(2,1,:) = -ep;
y(2,2,:) = -p1.*ep.*x1.^2;
y(2,3,:) = -p1.*ep.*x1.^4;

end

function y=obv_dxdx(x,p)

x1 = x(1,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
pp = x1+p2.*x1.^2+p3.*x1.^4;
pd = 1+2*p2.*x1+4*p3.*x1.^3;
pt = 2*p2+12*p3.*x1.^2;
ep = exp(pp);

y = zeros(2,2,2,numel(x1));
y(2,1,1,:)=-p1.*ep.*(pt+pd.^2);

end

function y=obv_dxdp(x,p)

x1 = x(1,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
pp = x1+p2.*x1.^2+p3.*x1.^4;
pd = 1+2*p2.*x1+4*p3.*x1.^3;
ep = exp(pp);

y = zeros(2,2,3,numel(x1));
y(2,1,1,:) = -ep.*pd;
y(2,1,2,:) = -ep.*p1.*x1.*(pd.*x1+2);
y(2,1,3,:) = -ep.*p1.*x1.^3.*(pd.*x1+4);

end

function y=obv_dpdp(x,p)

x1 = x(1,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
pp = x1+p2.*x1.^2+p3.*x1.^4;
ep = exp(pp);

y = zeros(2,3,3,numel(x1));
y(2,1,2,:) = -ep.*x1.^2;
y(2,1,3,:) = -ep.*x1.^4;
y(2,2,1,:) = -ep.*x1.^2;
y(2,2,2,:) = -p1.*ep.*x1.^4;
y(2,2,3,:) = -p1.*ep.*x1.^6;
y(2,3,1,:) = -ep.*x1.^4;
y(2,3,2,:) = -p1.*ep.*x1.^6;
y(2,3,3,:) = -p1.*ep.*x1.^8;

end

function [data, y] = obv_bc(prob, data, u) %#ok<INUSL>

T  = u(1);
x0 = u(2:3);
x1 = u(4:5);
p  = u(6:8); %#ok<NASGU>

y = [T-1; x0(1); x1(1)];

end

function [data, J] = obv_bc_du(prob, data, u) %#ok<INUSD,INUSL>

J = zeros(3,8);
J(1,1) = 1;
J(2,2) = 1;
J(3,4) = 1;

end

function [data, dJ] = obv_bc_dudu(prob, data, u) %#ok<INUSD,INUSL>

dJ = zeros(3,8,8);

end

function data = obj_init_data(fdata)

data.coll_seg = fdata.coll_seg;
data.ghan     = @ghan;
data.ghan_dx  = @ghan_dx;

data = coco_func_data(data);

end

function [prob, status, xtr] = obj_remesh(prob, data, chart, old_u, old_V) %#ok<INUSD>

[fdata, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
data.coll_seg = fdata.coll_seg;

xtr    = [];
prob   = coco_change_func(prob, data, 'uidx', uidx);
status = 'success';

end

function [data, y] = objhan(prob, data, u) %#ok<INUSL>

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T    = u(maps.T_idx);
x    = u(maps.xbp_idx);
p    = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
gcn = pr.ghan(xcn);
gcn = mesh.gka.*gcn;   % mesh here is related adaptive meshing

y = (0.5*T/maps.NTST)*mesh.gwt*gcn' + (p(1)^2+p(2)^2+p(3)^2)/10;

end

function [data, J] = objhan_du(prob, data, u) %#ok<INUSL>

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T = u(maps.T_idx);
x = u(maps.xbp_idx);
p = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
gcn = pr.ghan(xcn);
gcn = mesh.gka.*gcn;

gdxcn = pr.ghan_dx(xcn);
gdxcn = mesh.gdxka.*gdxcn;
gdxcn = sparse(maps.gdxrows, maps.gdxcols, gdxcn(:));

J_xbp = (0.5*T/maps.NTST)*mesh.gwt*gdxcn*maps.W;
J_T0  = 0;
J_T   = (0.5/maps.NTST)*mesh.gwt*gcn';
J_p   = [p(1) p(2) p(3)]/5;

J = [J_xbp J_T0 J_T J_p];

end

function data = adj_obj_init_data(fdata)

data.coll_seg  = fdata.coll_seg;
data.ghan      = @ghan;
data.ghan_dx   = @ghan_dx;
data.ghan_dxdx = @ghan_dxdx;

seg  = fdata.coll_seg;
maps = seg.maps;
int  = seg.int;

NCOL = int.NCOL;
NTST = maps.NTST;
xdim = int.dim;
pdim = maps.pdim;

rows = NCOL*NTST*kron(0:(xdim-1), ones(1,xdim));
opt.gdxdxrows1 = repmat(rows, [1 NCOL*NTST]) + ...
  kron(1:NCOL*NTST, ones(1,xdim^2));
cols = reshape(1:xdim*NCOL*NTST, [xdim NCOL*NTST]);
opt.gdxdxcols1 = repmat(cols, [xdim 1]);

% Derivative of (T/2N)*gxcn with respect to xbp:
step = 1+xdim*(0:NCOL-1); % NCOL
step = repmat(step(:), [1 xdim]) + repmat(0:xdim-1, [NCOL 1]); % xdim*NCOL
step = repmat(step(:), [1 xdim*(NCOL+1)]) + ...
  (xdim*NCOL*NTST+2+pdim)*repmat(0:xdim*(NCOL+1)-1, [xdim*NCOL 1]); % xdim^2*NCOL*(NCOL+1)
step = repmat(step(:), [1 NTST]) + ...
  (xdim*NCOL+xdim*(NCOL+1)*(xdim*NCOL*NTST+2+pdim))*...
  repmat(0:NTST-1, [xdim^2*NCOL*(NCOL+1) 1]); % xdim^2*NCOL*(NCOL+1)*NTST
opt.gdxdxcols2 = step(:);
opt.gdxdxrows2 = ones(xdim^2*NCOL*(NCOL+1)*NTST, 1);

step = 1:NCOL; % NCOL
step = repmat(step(:), [1 xdim]) + NTST*NCOL*repmat(0:xdim-1, [NCOL 1]); % xdim*NCOL
step = repmat(step(:), [1 xdim*(NCOL+1)]) + ...
  xdim*NTST*NCOL*repmat(0:xdim*(NCOL+1)-1, [xdim*NCOL 1]); % xdim^2*NCOL*(NCOL+1)
step = repmat(step(:), [1 NTST]) + (NCOL+xdim^2*NTST*NCOL*(NCOL+1))*...
  repmat(0:NTST-1, [xdim^2*NCOL*(NCOL+1) 1]); % xdim^2*NCOL*(NCOL+1)*NTST
opt.gdxdxidx = step(:);

% Derivative of (T/2N)*gxcn with respect to T:
opt.gdxdTrows = ones(xdim*NTST*NCOL, 1);
opt.gdxdTcols = (xdim*NCOL*NTST+2+pdim)*xdim*(NCOL+1)*NTST + ...
  (1:xdim*NTST*NCOL)';

% Derivative of (1/2N)*w*g' with respect to xbp:
opt.gdTdxcols = NTST*NCOL*xdim+2 + ...
  (xdim*NTST*NCOL+2+pdim)*(0:xdim*(NCOL+1)*NTST-1)';
opt.gdTdxrows = ones(xdim*(NCOL+1)*NTST, 1);

% Derivative of [p1/5, p2/5, p3/5] with respect to p:

opt.gdpdprows = ones(3,1);
opt.gdpdpcols = (NTST*xdim*NCOL+2+pdim)*(xdim*(NCOL+1)*NTST+2) + ...
  xdim*NTST*NCOL+2+[1 NTST*xdim*NCOL+2+pdim+2 2*(NTST*xdim*NCOL+2+pdim)+3]';

opt.dJrows = 1;
opt.dJcols = (xdim*NTST*NCOL+2+pdim)*(xdim*NTST*(NCOL+1)+2+pdim);

data.coll_opt = opt;

data = coco_func_data(data);

end

function [prob, status, xtr, ftr] = adj_obj_remesh(prob, data, chart, lb, Vlb)  %#ok<INUSL>

[fdata, axidx] = coco_get_adjt_data(prob, 'coll', 'data', 'axidx');
data = adj_obj_init_data(fdata);
fopt = fdata.coll_opt;

aidx = axidx([fopt.xcn_idx; fopt.T0_idx; fopt.T_idx; fopt.p_idx]);

xtr    = [];
ftr    = 1;
prob   = coco_change_adjt(prob, data, 'aidx', aidx, 'l0', lb, 'vecs', Vlb);
status = 'success';

end

function [data, y] = adj_objhan(prob, data, u) %#ok<INUSL>

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T = u(maps.T_idx);
x = u(maps.xbp_idx);
p = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
gcn = pr.ghan(xcn);
gcn = mesh.gka.*gcn;

gdxcn = pr.ghan_dx(xcn);
gdxcn = mesh.gdxka.*gdxcn;

J_xbp = (0.5*T/maps.NTST)*gdxcn(:)';
J_T0  = 0;
J_T   = (0.5/maps.NTST)*mesh.gwt*gcn';
J_p   = [p(1) p(2) p(3)]/5;

y = [J_xbp J_T0 J_T J_p];

end

function [data, J] = adj_objhan_du(prob, data, u) %#ok<INUSL>

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;
opt  = pr.coll_opt;

T = u(maps.T_idx);
x = u(maps.xbp_idx);

xcn = reshape(maps.W*x, maps.x_shp);

gdxdxcn = pr.ghan_dxdx(xcn);
gdxdxcn = mesh.gdxdxka.*gdxdxcn;
gdxdxcn = sparse(opt.gdxdxrows1, opt.gdxdxcols1, gdxdxcn(:))*maps.W;
J       = (0.5*T/maps.NTST)*sparse(opt.gdxdxrows2, opt.gdxdxcols2, ...
  gdxdxcn(opt.gdxdxidx), opt.dJrows, opt.dJcols);

gdxcn   = pr.ghan_dx(xcn);
gdxcn   = mesh.gdxka.*gdxcn;
J       = J + (0.5/maps.NTST)*sparse(opt.gdxdTrows, opt.gdxdTcols, ...
  gdxcn(:), opt.dJrows, opt.dJcols);

gdxcn   = mesh.gwt*sparse(maps.gdxrows, maps.gdxcols, gdxcn(:))*maps.W;
J       = J + (0.5/maps.NTST)*sparse(opt.gdTdxrows, opt.gdTdxcols, ...
  gdxcn(:), opt.dJrows, opt.dJcols);

J       = J + sparse(opt.gdpdprows, opt.gdpdpcols, ones(1,3)/5, ...
  opt.dJrows, opt.dJcols);

end

function y = ghan(x)

x1 = x(1,:);

y  = (x1-1).^2;

end

function J = ghan_dx(x)

x1 = x(1,:);

J  = zeros(1,2,numel(x1));
J(1,1,:) = 2*(x1-1);

end

function dJ = ghan_dxdx(x)

x1 = x(1,:);

dJ = zeros(1,2,2,numel(x1));
dJ(1,1,1,:) = 2;

end


function y = obv_ode(t,x,p1,p2,p3)

x1 = x(1);
x2 = x(2);
y = [x2; -p1*exp(x1+p2*x1.^2+p3*x1.^4)];

end