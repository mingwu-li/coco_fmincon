function [sol data] = ddaecoll_read_solution(oid, run, lab)
%COLL_READ_SOLUTION   Read 'coll' solution and toolbox data from disk.
%
% Extract data and chart structures associated with 'coll' toolbox instance
% identifier from solution file and construct solution structure.
%
% [SOL DATA] = COLL_READ_SOLUTION(OID, RUN, LAB)
%
% SOL  - Solution (struct).
% DATA - Toolbox data (struct).
% OID  - Object instance identifier (string).
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_read_solution.m 2839 2015-03-05 17:09:01Z fschild $

tbid         = coco_get_id(oid, 'ddaecoll');
[data chart] = coco_read_solution(tbid, run, lab);

sol.t  = chart.x(data.T0_idx)+data.tbpd(data.tbp_idx)*chart.x(data.T_idx);
xbp    = reshape(chart.x(data.xbp_idx), data.xbp_shp)';
sol.x  = xbp(data.tbp_idx,:);
sol.ta = chart.x(data.T0_idx)+data.tbpa*chart.x(data.T_idx);
ybp    = reshape(chart.x(data.ybp_idx), data.y_shp)';
sol.ya = ybp;
sol.p  = chart.x(data.p_idx);

if strcmp(data.ddaecoll.Apoints,'Uniform')
    ta    = data.tbpa(data.taL_idx);
    ya    = (ybp(data.taL_idx,:)+ybp(data.taR_idx,:))/2;
    sol.y = interp1(ta,ya,data.tbpd(data.tbp_idx),'pchip','extrap');
else
    sol.y = interp1(data.tbpa,ybp,data.tbpd(data.tbp_idx),'pchip','extrap');
end

end
