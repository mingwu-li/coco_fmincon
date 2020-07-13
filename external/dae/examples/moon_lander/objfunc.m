function [y, Dy] =  objfunc(u, prob, objfid)
    Fi = prob.efunc.funcs(strcmp(objfid,prob.efunc.identifyers));
    ui = u(Fi.x_idx);
    [~,y] = Fi.F(prob, Fi.data, ui);
    if nargout>1
       Dy = zeros(numel(u),1);
       if ~isempty(Fi.DFDX)
           [~, Dyi] = Fi.DFDX(prob, Fi.data, ui);
       elseif nargout(Fi.F)==3
           [~, ~, Dyi] = Fi.F(prob, Fi.data, ui);
       else
           [~, Dyi] = coco_ezDFDX('f(o,d,x)', prob, Fi.data, Fi.F, ui);
       end
       Dy(Fi.x_idx) = Dyi(:);
    end
end