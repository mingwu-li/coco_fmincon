function [x,y] = opt_read_sol(u, prob, fid)
% This function returns the design variables specific to a function with
% identifyer fid, and the evaluation of this function at such variables.
%
% Here the input u is the vector of all design variables, the prob is a
% problem instance of coco, and fid is the identifyer of the function.
% In the output, x represents the design variables specific to the function
% and y gives the evaluation of the function.

%%
% find the function from prob based on its identifyer
Fi = prob.efunc.funcs(strcmp(fid, prob.efunc.identifyers));
% extract x from u
x  = u(Fi.x_idx);
% evaluation of Fi at x
[~,y] = Fi.F(prob, Fi.data, x); 

end