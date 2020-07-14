function [y, Dy] =  objfunc(u, prob, objfid)
% This function constructs the objective to be minimized from the prob of
% coco. 
%
% Here the input u is a vector of design variables, the prob is a
% problem instance of coco, objfid is a string for the function identifyer
% of the objective function(al).
% In the output, y is a scalar for the objective function(al), and Dy is
% its gradient.

%%
% find objetive function(al) based on its identifyer
Fi = prob.efunc.funcs(strcmp(objfid,prob.efunc.identifyers));
ui = u(Fi.x_idx);
[~,y] = Fi.F(prob, Fi.data, ui);
assert(numel(y)==1, 'Objetive is not scalar-valued function');
% gradient
if nargout>1
   Dy = zeros(numel(u),1);
   if ~isempty(Fi.DFDX)      % with nonempty DFDX
       [~, Dyi] = Fi.DFDX(prob, Fi.data, ui);
   elseif nargout(Fi.F)==3   % with [o,y,J] as output to Fi
       [~, ~, Dyi] = Fi.F(prob, Fi.data, ui);
   else                      % numerical differentiation
       [~, Dyi] = coco_ezDFDX('f(o,d,x)', prob, Fi.data, Fi.F, ui);
   end
   Dy(Fi.x_idx) = Dyi(:);
end
    
    
end