function sol = ddaecoll_init_sol(data, t0, x0, y0, p0)
%COLL_INIT_SOL   Build initial solution guess.
%
% Use sampled trajectory on temporal mesh to construct initial solution
% guess for 'coll' toolbox.
%
% SOL = COLL_INIT_SOL(DATA, T0, X0, P0)
%
% DATA - Toolbox data structure.
% T0   - Array of temporal mesh points.
% X0   - Array of state vectors at mesh points.
% P0   - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_init_sol.m 2839 2015-03-05 17:09:01Z fschild $

t0 = t0(:);
T0 = t0(1);
T  = t0(end)-t0(1); % Interval length
t0 = (t0-t0(1))/T;  % Rescaling - only if T0<>0!
x0 = interp1(t0, x0, data.tbpd)'; % Collection of basepoint values
if ~isempty(y0)
y0 = interp1(t0, y0, data.tbpa)';
end
sol.u = [x0(:); y0(:); T0; T; p0(:)];

end
