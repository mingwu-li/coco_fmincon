function [c, y, Dc, Dy] =  nonlincons(u, prob)
% This function constructs inequality and equality constraints of
% an optimization problem from the prob of coco. Equality constraints are
% defined as zero functions. So all zero functions are used to construct
% corresponding equality constraints. Inequality constraints are defined as
% inequality monitor functions. So all inequality functions are used to
% construct corresponding inequality constraints.
%
% Here the input u is a vector of design variables, and the prob is a
% problem instance of coco.
% In the output, c respresents inequality constraints (c<=0) and y gives
% equality constraints (y=0). Dc and Dy is the Jacobian of c and y with
% respect to u respectively.
%
% For gradient-based optimization methods, the user may provide such
% Jacobian and specify proper option for gradient in the setting of
% fmincon.

%%
% inequality
c = [];
for i = prob.efunc.embedded
    Fi = prob.efunc.funcs(i);
    if strcmp(Fi.type, 'inequality') % inequality functions
        ui = u(Fi.x_idx);
        [~, ci] = Fi.F(prob, Fi.data, ui); 
        c = [c; ci(:)];
    end
end

% equality
y = cell(1,length(prob.efunc.zero));
for i=prob.efunc.zero                % zero functions
    Fi = prob.efunc.funcs(i);
    ui = u(Fi.x_idx);
    [~,y{i}] = Fi.F(prob, Fi.data, ui);
end
y = cell2mat(y(:));

% Gradient of the constraints:
if nargout > 2
    % inequality
    Dc = zeros(numel(c),numel(u));
    cdim = 0;
    for i = prob.efunc.embedded
        Fi  = prob.efunc.funcs(i);
        if strcmp(Fi.type, 'inequality')
            ui  = u(Fi.x_idx);
            Fi_dim = numel(Fi.f_idx);
            Fi_idx = cdim+1:cdim+Fi_dim;
            if ~isempty(Fi.DFDX)     % with nonempty DFDX
                [~, Dc(Fi_idx,Fi.x_idx)] = Fi.DFDX(prob, Fi.data, ui);
            elseif nargout(Fi.F)==3  % with [o,y,J] as output to Fi
                [~, ~, Dc(Fi_idx,Fi.x_idx)] = Fi.F(prob, Fi.data, ui);
            else                     % numerical differentiation
                [~, Dc(Fi_idx,Fi.x_idx)] = coco_ezDFDX('f(o,d,x)', prob, Fi.data, Fi.F, ui);
            end
            cdim = cdim + Fi_dim;
        end
    end  
    Dc = Dc'; % fmincon asks each column of Dy corresponds to the gradient of one constraint

    % equality
    Dy = zeros(numel(y),numel(u));
    Fdim = 0;
    for i = prob.efunc.zero
        Fi  = prob.efunc.funcs(i);
        ui  = u(Fi.x_idx);
        Fi_dim = numel(Fi.f_idx);
        Fi_idx = Fdim+1:Fdim+Fi_dim;
        if ~isempty(Fi.DFDX)     % with nonempty DFDX
            [~, Dy(Fi_idx,Fi.x_idx)] = Fi.DFDX(prob, Fi.data, ui);
        elseif nargout(Fi.F)==3  % with [o,y,J] as output to Fi
            [~, ~, Dy(Fi_idx,Fi.x_idx)] = Fi.F(prob, Fi.data, ui);
        else                     % numerical differentiation
            [~, Dy(Fi_idx,Fi.x_idx)] = coco_ezDFDX('f(o,d,x)', prob, Fi.data, Fi.F, ui);
        end
        Fdim = Fdim + Fi_dim;
    end  
    Dy = Dy'; % fmincon asks each column of Dy corresponds to the gradient of one constraint
end

end


