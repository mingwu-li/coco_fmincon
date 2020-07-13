function prob = alg_dae_construct_seg(prob, tbid, data)
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

assert(numel(data.uidx)==data.ddaecoll_seg.T_idx+numel(data.ddaecoll_seg.p_idx), ...
    'size of uidx do not match');
if isempty(data.ineq) || ~strcmpi(data.ineq, 'inequality')
    prob = coco_add_func(prob, tbid, @alg_dae_F, @alg_dae_DFDU, data, 'zero', ...
      'uidx', data.uidx);
else  
    seg  = data.ddaecoll_seg;
    dim  = seg.dim;    % State-space dimension
    ydim = seg.ydim;
    pdim = seg.pdim;      % Number of problem parameters
    gdim = numel(data.ghan(0, zeros(dim,1), zeros(ydim,1), zeros(pdim,1))); % not necessary here given the output of g is of the same dimension as pdim
    ineqid  = cell(1, gdim*data.ycnnum);
    ineqmid = coco_get_id(tbid, 'ineq');
    for k=1:gdim*data.ycnnum
        ineqid{k} = coco_get_id(ineqmid, num2str(k));
    end
    prob = coco_add_func(prob, tbid, @alg_dae_F, @alg_dae_DFDU, data, 'inequality', ...
      ineqid, 'uidx', data.uidx);
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

function [data y] = alg_dae_F(prob, data, u)
%COLL_F   Collocation zero problem.
%
% Expects vectorized encoding of vector field as function of two arguments.

seg = data.ddaecoll_seg;
x   = u(seg.xbp_idx); % Extract basepoint values
y   = u(seg.ybp_idx);
T0  = u(seg.T0_idx);
T   = u(seg.T_idx);   % Extract interval length
p   = u(seg.p_idx);   % Extract problem parameters

xx  = reshape(seg.Wad*x, seg.x_shp); % Values at collocation nodes
pp  = repmat(p, seg.p_rep);
yy  = reshape(seg.Waa*y, seg.y_shp); 
tcn = seg.Ma*T+T0;
alg = data.ghan(tcn',xx, yy, pp);

y = alg(:);

end

function [data J] = alg_dae_DFDU(prob, data, u)
%COLL_DFDU   Linearization of collocation zero problem.
%
% Expects vectorized encoding of Jacobians of vector field with respect to
% its arguments. When dfdxhan and/or dfdphan are empty, approximate
% Jacobians are obtained using numerical differentiation.

seg = data.ddaecoll_seg;
x   = u(seg.xbp_idx); % Extract basepoint values
y   = u(seg.ybp_idx);
T0  = u(seg.T0_idx);
T   = u(seg.T_idx);   % Extract interval length
p   = u(seg.p_idx);   % Extract problem parameters

xx  = reshape(seg.Wad*x, seg.x_shp); % Values at collocation nodes
pp  = repmat(p, seg.p_rep);
yy  = reshape(seg.Waa*y, seg.y_shp); 
tcn = seg.Ma*T+T0;

ydim = seg.ydim;
dim  = seg.dim;

if isempty(data.dgdxhan)
  f  = @(x,p) data.ghan(p(ydim+1,:), x, p(1:ydim,:), p(ydim+2:end,:));
  ytp = [yy; tcn'; pp];
  dxode = coco_ezDFDX('f(x,p)v', f, xx, ytp);
else
  dxode = data.dgdxhan(tcn', xx, yy, pp);
end
dxode = sparse(data.gdxrows, data.gdxcols, dxode(:));
dxode = dxode*seg.Wad; % W.r.t. basepoint values

if isempty(data.dgdyhan)
  f  = @(x,p) data.ghan(p(dim+1,:), p(1:dim,:), x, p(dim+2:end,:));
  xtp = [xx; tcn'; pp];
  dyode = coco_ezDFDX('f(x,p)v', f, yy, xtp);    
else
  dyode = data.dgdyhan(tcn', xx, yy, pp);
end
dyode = sparse(data.gdyrows, data.gdycols, dyode(:));
dyode = dyode*seg.Waa;

if isempty(data.dgdthan)
  f  = @(x,p) data.ghan(x, p(1:dim,:), p(dim+1:dim+ydim,:), p(dim+ydim+1:end,:));
  xyp = [xx; yy; pp];
  dtode = coco_ezDFDX('f(x,p)v', f, tcn', xyp);
else
  dtode = data.dgdthan(tcn', xx, yy, pp);
end
dtode  = sparse(data.gdtrows, data.gdtcols, dtode(:));
dT0ode = dtode*ones(data.ycnnum,1);
dTode  = dtode*seg.Ma; 

if isempty(data.dgdphan)
  f  = @(x,p) data.ghan(p(dim+ydim+1,:), p(1:dim,:), p(dim+1:dim+ydim,:), x);
  xyt = [xx; yy; tcn'];
  dpode = coco_ezDFDX('f(x,p)v', f, pp, xyt);
else
  dpode = data.dgdphan(tcn', xx, yy, pp);
end
dpode = sparse(data.gdprows, data.gdpcols, dpode(:)); % W.r.t. problem parameters

J = [dxode dyode dT0ode dTode dpode]; 
% [data, DadF] = coco_ezDFDX('f(o,d,x)',  prob, data, @alg_dae_F, u);

end
