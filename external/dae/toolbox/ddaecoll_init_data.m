function data = ddaecoll_init_data(data, x0, y0, p0)
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

NTST = data.ddaecoll.NTST; % Number of mesh intervals
NCOL = data.ddaecoll.NCOL; % Degree of polynomial interpolants
dim  = size(x0, 2);    % State-space dimension
ydim = size(y0,2);
pdim = numel(p0);      % Number of problem parameters

data.dim  = dim;       % State-space dimension
data.ydim = ydim;      % Dimension of y variables
data.pdim = pdim;      % Number of problem parameters

bpnum  = NCOL+1;            % Number of basepoints per interval
bpdim  = dim*(NCOL+1);      % Number of basepoint values per interval
xbpnum = (NCOL+1)*NTST;     % Number of basepoints
xbpdim = dim*(NCOL+1)*NTST; % Number of basepoint values
ybpdim = ydim*NCOL*NTST;    % Number of basepoint values for y
cndim  = dim*NCOL;          % Number of collocation node values per interval
xcnnum = NCOL*NTST;         % Number of collocation nodes
xcndim = dim*NCOL*NTST;     % Number of collocation conditions
ycndim = ybpdim;            % Number of collocation conditions for y 
cntnum = NTST-1;            % Number of internal boundaries
cntdim = dim*(NTST-1);      % Number of continuity conditions

data.xbpnum  = xbpnum;
data.xcnnum  = xcnnum;
data.xbp_idx = (1:xbpdim)'; % Index array for basepoint values
data.ybp_idx = xbpdim+(1:ybpdim)';
data.T0_idx  = xbpdim+ybpdim+1;
data.T_idx   = xbpdim+ybpdim+2;    % Index for interval length
data.p_idx   = xbpdim+ybpdim+2+(1:pdim)'; % Index array for problem parameters
data.tbp_idx = setdiff(1:xbpnum, 1+bpnum*(1:cntnum))'; % Index array without duplication of internal boundaries
data.taL_idx = setdiff(1:xcnnum, 1+NCOL*(1:cntnum))';  % Index array without duplication for uniform algebraic basepoints.
data.taR_idx = setdiff(1:xcnnum, NCOL*(1:cntnum))';    % Index array without duplication for uniform algebraic basepoints.
data.x_shp   = [dim xcnnum]; % Shape for vectorization
data.xbp_shp = [dim xbpnum]; % Shape for vectorization
data.y_shp   = [ydim xcnnum];
data.p_rep   = [1 xcnnum];   % Shape for vectorization


switch data.ddaecoll.Dpoints
    case 'Uniform'  % default setting
        tsd = linspace(-1, 1, NCOL+1)';
    case 'Gauss_Radau'
        tsd = [legsrd(NCOL);1];
    case 'Flipped_Gauss_Radau'
        tsd = legsrd(NCOL); 
        tsd = [-1; -tsd(end:-1:1)];
    otherwise
        disp('Please select Dpoints from {Uniform,Gauss_Radau,Flipped_Gauss_Radau}');
end

switch data.ddaecoll.Apoints
    case 'Uniform'  % default setting
        tsa = linspace(-1,1,NCOL)';
    case 'Gauss'  
        [tsa,~] = coll_nodes(NCOL);
    case 'Gauss_Radau'
        tsa = legsrd(NCOL);
    case 'Flipped_Gauss_Radau'
        tsa = legsrd(NCOL);
        tsa = -tsa(end:-1:1);
    otherwise
        disp('Please select Apoints from {Uniform,Gauss,Gauss_Radau,Flipped_Gauss_Radau}');
end

switch data.ddaecoll.Dnodes
    case 'Gauss'  % default setting
        [tzd, wts] = coll_nodes(NCOL);
    case 'Gauss_Radau'
        [tzd, wts] = legsrd(NCOL); 
    case 'Flipped_Gauss_Radau'
        [tzd, wts] = legsrd(NCOL); 
        tzd = -tzd(end:-1:1);
        wts  = wts(end:-1:1);
    otherwise
        disp('Please select Dnodes from {Uniform,Gauss,Gauss_Radau,Flipped_Gauss_Radau}');
end

switch data.ddaecoll.Anodes
    case 'Uniform'
        tza = linspace(-1,1,NCOL)';
    case 'Gauss'  % default setting
        [tza, ~] = coll_nodes(NCOL);
    case 'Gauss_Radau'
        tza = legsrd(NCOL); 
    case 'Flipped_Gauss_Radau'
        tza = legsrd(NCOL);
        tza = -tza(end:-1:1);
    otherwise
        disp('Please select from {Uniform,Gauss,Gauss_Radau,Flipped_Gauss_Radau}');
end

data.tsa    = tsa;
data.tsd    = tsd;
data.tza    = tza;
data.tzd    = tzd;

tbpd         = kron((0:NTST-1)',ones(NCOL+1,1))+repmat(0.5*(tsd+1),[NTST,1]);
data.tbpd    = tbpd/NTST;
tbpa         = kron((0:NTST-1)',ones(NCOL,1))+repmat(0.5*(tsa+1),[NTST,1]);
data.tbpa    = tbpa/NTST;

Md          = kron((0:NTST-1)',ones(NCOL,1))+repmat(0.5*(tzd+1),[NTST,1]);
data.Md     = Md/NTST;
Ma          = kron((0:NTST-1)',ones(NCOL,1))+repmat(0.5*(tza+1),[NTST,1]);
data.Ma     = Ma/NTST;

data.x0_idx = (1:dim)';            % Index array for trajectory end point at t=0
data.x1_idx = xbpdim-dim+(1:dim)'; % Index array for trajectory end point at t=1

data.wts    = wts;
wts         = repmat(wts, [dim NTST]);
data.wts1   = wts(1,:);                          
data.wts2   = spdiags(wts(:), 0, xcndim, xcndim);

Ldd         = coll_L(tsd, tzd);
Lddp        = coll_Lp(tsd, tzd);
rows        = reshape(1:xcndim, [cndim NTST]);
rows        = repmat(rows, [bpdim 1]);
cols        = repmat(1:xbpdim, [cndim 1]);
Wdd         = repmat(kron(Ldd, eye(dim)), [1 NTST]);  % Interpolation matrix
Wddp        = repmat(kron(Lddp, eye(dim)), [1 NTST]); % Interpolation matrix
data.Wdd    = sparse(rows, cols, Wdd);
data.Wddp   = sparse(rows, cols, Wddp);

Lad         = coll_L(tsd, tza);
Wad         = repmat(kron(Lad, eye(dim)), [1 NTST]);  % Interpolation matrix
data.Wad    = sparse(rows, cols, Wad);

Lda         = coll_L(tsa, tzd);
yscndim     = ydim*NCOL; % Number of collocation nodes per subinterval
rows        = reshape(1:ycndim, [yscndim NTST]);
rows        = repmat(rows, [yscndim 1]);
cols        = repmat(1:ybpdim, [yscndim 1]);
Wda         = repmat(kron(Lda, eye(ydim)), [1 NTST]); % Interpolation matrix
data.Wda    = sparse(rows, cols, Wda); 

Laa         = coll_L(tsa, tza);
Waa         = repmat(kron(Laa, eye(ydim)), [1 NTST]); % Interpolation matrix
data.Waa    = sparse(rows, cols, Waa);

rows        = reshape(1:xcndim, [cndim NTST]);
rows        = repmat(rows, [cndim 1]);
cols        = repmat(1:xcndim, [cndim 1]);
Wyy         = repmat(kron(Laa, eye(dim)), [1 NTST]); % Interpolation matrix
data.Wyy    = sparse(rows, cols, Wyy);               % For DDEs specifically

data.dxrows = repmat(reshape(1:xcndim, [dim xcnnum]), [dim 1]);  % Index array for vectorization
data.dxcols = repmat(1:xcndim, [dim 1]);                         % Index array for vectorization
data.dprows = repmat(reshape(1:xcndim, [dim xcnnum]), [pdim 1]); % Index array for vectorization
data.dpcols = repmat(1:pdim, [dim xcnnum]);                      % Index array for vectorization
data.dyrows = repmat(reshape(1:xcndim, [dim xcnnum]), [ydim 1]);  % Index array for vectorization
data.dycols = repmat(1:ycndim, [dim 1]);                          % Index array for vectorization
data.dtrows = reshape(1:xcndim, [dim xcnnum]);                    % Index array for vectorization
data.dtcols = repmat(1:xcnnum, [dim 1]);                          % Index array for vectorization

temp        = reshape(1:xbpdim, [bpdim NTST]);
Qrows       = [1:cntdim 1:cntdim];
Qcols       = [temp(1:dim, 2:end) temp(cndim+1:end, 1:end-1)];
Qvals       = [ones(cntdim,1) -ones(cntdim,1)];
data.Q      = sparse(Qrows, Qcols, Qvals, cntdim, xbpdim); % Jacobian of continuity conditions
data.dyTpcnt= sparse(cntdim, ycndim+2+pdim);

end

function [nds wts] = coll_nodes(m)
%COLL_NODES   Compute collocation nodes and integration weights.
%
% Uses eigenvalues and eigenvectors of Jacobi matrix.
%
% [NDS WTS] = COLL_NODES(M)
%
% NDS - Collocation nodes.
% WTS - Quadrature weights.
% M   - Polynomial degree.

n = (1:m-1)';
g = n.*sqrt(1./(4*n.^2-1));
J = -diag(g,1)-diag(g,-1);

[w x] = eig(J);
nds   = diag(x);
wts   = 2*w(1,:).^2;

end

function A = coll_L(ts, tz)
%COLL_L   Evaluation of Lagrange polynomials.
%
% Use high-dimensional arrays for vectorized evaluation.
%
% A = COLL_L(TS, TZ)
% 
% A  - Array of interpolated values.
% TS - Array of basepoints.
% TZ - Array of interpolation points.

q = numel(ts);
p = numel(tz);

zi = repmat(reshape(tz, [p 1 1]), [1 q q]);
sj = repmat(reshape(ts, [1 q 1]), [p 1 q]);
sk = repmat(reshape(ts, [1 1 q]), [p q 1]);

t1 = zi-sk;
t2 = sj-sk;
idx = find(abs(t2)<=eps);
t1(idx) = 1;
t2(idx) = 1;

A = prod(t1./t2, 3);

end

function A = coll_Lp(ts, tz)
%COLL_LP   Evaluation of derivative of Lagrange polynomials.
%
% Use high-dimensional arrays for vectorized evaluation.
%
% A = COLL_LP(TS, TZ)
% 
% A  - Array of interpolated values
% TS - Array of basepoints
% TZ - Array of interpolation points

q = numel(ts);
p = numel(tz);

zi = repmat(reshape(tz, [p 1 1 1]), [1 q q q]);
sj = repmat(reshape(ts, [1 q 1 1]), [p 1 q q]);
sk = repmat(reshape(ts, [1 1 q 1]), [p q 1 q]);
sl = repmat(reshape(ts, [1 1 1 q]), [p q q 1]);

t3 = sj(:,:,:,1)-sk(:,:,:,1);
t4 = zi-sl;
t5 = sj-sl;

idx1 = find(abs(t5)<=eps);
idx2 = find(abs(t3)<=eps);
idx3 = find(abs(sk-sl)<=eps);
t5(union(idx1, idx3)) = 1;
t4(union(idx1, idx3)) = 1;
t3(idx2) = 1;
t3       = 1.0./t3;
t3(idx2) = 0;

A = sum(t3.*prod(t4./t5, 4), 3);

end



 function [varargout]=legsrd(n)

 % x=legsrd(n) returns n Legendre-Gauss-Radau points with x(1)=-1.
%  [x,w]= legsrd(n) returns n Legendre-Gauss-Radau points and weights
%  Eigenmethod is used for computing nodes.
%  Last modified on August 30, 2011

 j=[0:n-2];                  % indices    
 A=diag(1./((2*j+1).*(2*j+3)));   % Main diagonal
 j=[1:n-2];
 A=A+diag(sqrt(j.*(j+1))./(2*j+1),1) ...
  +diag(sqrt(j.*(j+1))./(2*j+1),-1);     %  Create Jacobi matrix
 x= sort(eig(sparse(A)));              %  Compute eigenvalues
 x=[-1;x];
     varargout{1}=x;
  if nargout==1, return; end;
y=legendreP(n-1,x);
varargout{2}= (1-x)./(n^2*y.^2);  % return the weights by (3.179)

 end
