function prob = ddaecoll_construct_seg(prob, tbid, data, sol)
%COLL_CONSTRUCT_SEG   Append an instance of 'coll' to problem.
%
% Add collocation and continuity conditions, monitor functions that
% evaluate to the problem parameters, and corresponding inactive
% continuation parameters.
%
% PROB = COLL_CONSTRUCT_SEG(PROB, TBID, DATA, SOL)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% SOL  - Initial solution guess.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_construct_seg.m 2839 2015-03-05 17:09:01Z fschild $

prob = coco_add_func(prob, tbid, @ddaecoll_F, @ddaecoll_DFDU, data, 'zero', ...
  'u0', sol.u);
if ~isempty(data.pnames) % Optional monitor functions
  fid  = coco_get_id(tbid, 'pars');
  uidx = coco_get_func_data(prob, tbid, 'uidx'); % Context-dependent index array
  prob = coco_add_pars(prob, fid, uidx(data.p_idx), data.pnames);
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

function [data y] = ddaecoll_F(prob, data, u)
%COLL_F   Collocation zero problem.
%
% Expects vectorized encoding of vector field as function of two arguments.

x  = u(data.xbp_idx); % Extract basepoint values
y  = u(data.ybp_idx);
T0 = u(data.T0_idx);
T  = u(data.T_idx);   % Extract interval length
p  = u(data.p_idx);   % Extract problem parameters

xx = reshape(data.Wdd*x, data.x_shp); % Values at collocation nodes
pp = repmat(p, data.p_rep);
yy = reshape(data.Wda*y, data.y_shp); 

tcn = data.Md*T+T0;

ode = data.fhan(tcn',xx, yy, pp);
ode = (0.5*T/data.ddaecoll.NTST)*ode(:)-data.Wddp*x;     % Collocation conditions
cnt = data.Q*x;                                % Continuity conditions

y = [ode; cnt];

end

function [data J] = ddaecoll_DFDU(prob, data, u)
%COLL_DFDU   Linearization of collocation zero problem.
%
% Expects vectorized encoding of Jacobians of vector field with respect to
% its arguments. When dfdxhan and/or dfdphan are empty, approximate
% Jacobians are obtained using numerical differentiation.

x  = u(data.xbp_idx); % Extract basepoint values
y  = u(data.ybp_idx);
T0 = u(data.T0_idx);
T  = u(data.T_idx);   % Extract interval length
p  = u(data.p_idx);   % Extract problem parameters

xx = reshape(data.Wdd*x, data.x_shp); % Values at collocation nodes
pp = repmat(p, data.p_rep);
yy = reshape(data.Wda*y, data.y_shp); 

tcn  = data.Md*T+T0;
ydim = data.ydim;
dim  = data.dim;

if isempty(data.dfdxhan)
  f  = @(x,p) data.fhan(p(ydim+1,:), x, p(1:ydim,:), p(ydim+2:end,:));
  ytp = [yy; tcn'; pp];
  dxode = coco_ezDFDX('f(x,p)v', f, xx, ytp);
else
  dxode = data.dfdxhan(tcn', xx, yy, pp);
end
dxode = sparse(data.dxrows, data.dxcols, dxode(:));
dxode = (0.5*T/data.ddaecoll.NTST)*dxode*data.Wdd-data.Wddp; % W.r.t. basepoint values

if isempty(data.dfdyhan)
  f  = @(x,p) data.fhan(p(dim+1,:), p(1:dim,:), x, p(dim+2:end,:));
  xtp = [xx; tcn'; pp];
  dyode = coco_ezDFDX('f(x,p)v', f, yy, xtp);    
else
  dyode = data.dfdyhan(tcn', xx, yy, pp);
end
dyode = (0.5*T/data.ddaecoll.NTST)*sparse(data.dyrows, data.dycols, dyode(:));
dyode = dyode*data.Wda;

if isempty(data.dfdthan)
  f  = @(x,p) data.fhan(x, p(1:dim,:), p(dim+1:dim+ydim,:), p(dim+ydim+1:end,:));
  xyp = [xx; yy; pp];
  dtode = coco_ezDFDX('f(x,p)v', f, tcn', xyp);
else
  dtode = data.dfdthan(tcn', xx, yy, pp);
end
dtode  = sparse(data.dtrows, data.dtcols, dtode(:));
dt0ode = (0.5*T/data.ddaecoll.NTST)*dtode;
dT0ode = dt0ode*ones(data.xcnnum,1);

dTode = data.fhan(tcn', xx, yy, pp);
dTode = (0.5/data.ddaecoll.NTST)*dTode(:)+dt0ode*data.Md; 

if isempty(data.dfdphan)
  f  = @(x,p) data.fhan(p(dim+ydim+1,:), p(1:dim,:), p(dim+1:dim+ydim,:), x);
  xyt = [xx; yy; tcn'];
  dpode = coco_ezDFDX('f(x,p)v', f, pp, xyt);
else
  dpode = data.dfdphan(tcn', xx, yy, pp);
end
dpode = sparse(data.dprows, data.dpcols, dpode(:));
dpode = (0.5*T/data.ddaecoll.NTST)*dpode; % W.r.t. problem parameters

J = [dxode dyode dT0ode dTode dpode; data.Q data.dyTpcnt]; % data.dyTpcnt = [0...0]
% [data, DadF] = coco_ezDFDX('f(o,d,x)',  prob, data, @ddaecoll_F, u);

end
