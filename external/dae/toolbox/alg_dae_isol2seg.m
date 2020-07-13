function prob = alg_dae_isol2seg(prob, oid, varargin)
%COLL_ISOL2SEG   Append 'coll' instance constructed from initial data.
%
% Parse input sequence to construct toolbox data and initial solution guess
% and use this to construct an instance of 'coll'.
%
% PROB     = COLL_ISOL2SEG(PROB, OID, VARARGIN)
% VARARGIN = { @F [(@DFDX | '[]') (@DFDy | '[]') [(@DFDP | '[]')]] T0 X0 [PNAMES] P0 }
%
% PROB   - Continuation problem structure.
% OID    - Object instance identifier (string).
% @F     - Function handle to vector field.
% @DFDX  - Optional function handle to Jacobian w.r.t. problem variables.
% @DFDP  - Optional function handle to Jacobian w.r.t. problem parameters.
% T0     - Array of temporal mesh points.
% X0     - Array of state vectors at mesh points (rectangular matrix with
%          size(x0,1) == numel(t0).
% PNAMES - Optional string label or cell array of string labels for
%          continuation parameters tracking problem parameters.
% P0     - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_isol2seg.m 2839 2015-03-05 17:09:01Z fschild $


tbid = coco_get_id(oid, 'alg_dae'); % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing

if ischar(str.peek)
    fbid = coco_get_id(str.get, 'ddaecoll');
else
    fbid = coco_get_id(oid, 'ddaecoll'); % Create toolbox instance identifier
end
[fdata, uidx]     = coco_get_func_data(prob, fbid, 'data', 'uidx');
data.ddaecoll_seg = fdata;
data.uidx         = uidx;

data.ghan = str.get;
data.dgdthan  = [];
data.dgdxhan  = [];
data.dgdyhan  = [];
data.dgdphan  = [];
data.dgdtdthan  = [];
data.dgdtdxhan  = [];
data.dgdtdyhan  = [];
data.dgdtdphan  = [];
data.dgdxdxhan  = [];
data.dgdxdyhan  = [];
data.dgdxdphan  = [];
data.dgdydyhan  = [];
data.dgdydphan  = [];
data.dgdpdphan  = [];
data.ineq       = [];

if is_empty_or_func(str.peek)
  data.dgdthan = str.get;
  if is_empty_or_func(str.peek)
    data.dgdxhan = str.get;
    if is_empty_or_func(str.peek)
        data.dgdyhan = str.get; 
        if is_empty_or_func(str.peek)
            data.dgdphan = str.get;
            if is_empty_or_func(str.peek)
                data.dgdtdthan = str.get;
                if is_empty_or_func(str.peek)
                    data.dgdtdxhan = str.get;
                    if is_empty_or_func(str.peek)
                        data.dgdtdyhan = str.get;
                        if is_empty_or_func(str.peek)
                            data.dgdtdphan = str.get;
                            if is_empty_or_func(str.peek)
                                data.dgdxdxhan = str.get;
                                if is_empty_or_func(str.peek)
                                    data.dgdxdyhan = str.get;
                                    if is_empty_or_func(str.peek)
                                        data.dgdxdphan = str.get;
                                        if is_empty_or_func(str.peek)
                                            data.dgdydyhan = str.get;
                                            if is_empty_or_func(str.peek)
                                                data.dgdydphan = str.get;
                                                if is_empty_or_func(str.peek)
                                                    data.dgdpdphan = str.get;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
  end
end

if is_empty_or_string(str.peek)
    ineq = str.get;
    if strcmpi('inequality',ineq)
        data.ineq = 'inequality';
    end
end

data = alg_dae_init_data(data);        % Build toolbox data
prob = alg_dae_construct_seg(prob, tbid, data); % Append continuation problem

end

function flag = is_empty_or_func(x)
%IS_EMPTY_OR_FUNC   Check if input is empty or contains a function handle.
flag = isempty(x) || isa(x, 'function_handle');
end


function flag = is_empty_or_string(x)
%IS_EMPTY_OR_FUNC   Check if input is empty or contains a function handle.
flag = isempty(x) || ischar(x);
end
