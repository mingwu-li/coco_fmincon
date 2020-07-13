function data = alg_dae_init_data(data)
%COLL_INIT_DATA   Initialize toolbox data for an instance of 'coll'.
%
% Populate remaining fields of the toolbox data structure used by 'coll'
% function objects.
%
% DATA = COLL_INIT_DATA(DATA, X0, P0)
%
% DATA - Toolbox data structure.
% X0   - Initial solution guess for discretized trajectory.
% P0   - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_init_data.m 2839 2015-03-05 17:09:01Z fschild $

seg  = data.ddaecoll_seg;
NTST = seg.ddaecoll.NTST; % Number of mesh intervals
NCOL = seg.ddaecoll.NCOL; % Degree of polynomial interpolants
dim  = seg.dim;    % State-space dimension
ydim = seg.ydim;
pdim = seg.pdim;      % Number of problem parameters

ycndim = ydim*NCOL*NTST;    % Number of basepoint values for y
ycnnum = NCOL*NTST;         % Number of collocation nodes
xcndim = dim*NCOL*NTST;     % Number of collocation conditions

data.ycnnum = ycnnum;

data.gdxrows = repmat(reshape(1:ycndim, [ydim ycnnum]), [dim 1]);   % Index array for vectorization
data.gdxcols = repmat(1:xcndim, [ydim 1]);                          % Index array for vectorization
data.gdprows = repmat(reshape(1:ycndim, [ydim ycnnum]), [pdim 1]);  % Index array for vectorization
data.gdpcols = repmat(1:pdim, [ydim ycnnum]);                       % Index array for vectorization
data.gdyrows = repmat(reshape(1:ycndim, [ydim ycnnum]), [ydim 1]);  % Index array for vectorization
data.gdycols = repmat(1:ycndim, [ydim 1]);                          % Index array for vectorization
data.gdtrows = reshape(1:ycndim, [ydim ycnnum]);                    % Index array for vectorization
data.gdtcols = repmat(1:ycnnum, [ydim 1]);                          % Index array for vectorization

end


