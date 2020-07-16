# A wrapper of COCO for fmincon in MATLAB

## COCO
Continuation Core (COCO) is a MATLAB-based open-source package for computational nonlinear analysis of dynamical systems. It uses atlas algorithms to cover the solution manifolds of nonlinear systems. A unique feature of COCO is the embedded construction philosphy of nonlinear systems, where a large problem is assembed by small subproblems with weak couplings. This object-oriented construction paradigm enables us to build a composite problem from building blocks. COCO provides constructors for collocation problems (coll), boundary value problems (bvp), periodic orbit problems (po) and hybrid system periodic orbit problems (hspo). Further, COCO provides constructors for the adjoint operators of these problems. Such adjoints can be used to formulate the first-order necessary conditions to constrained optimization problems. Please refer [1-4] for more info of COCO.


## fmincon
**fmincon** is a MATLAB function for finding the minimum of *finite-dimensional* constrained optimization problem. It is a nonlinear programming solver. More details about **fmincon** can be found at MATLAB documentation [5]. Note that **fmincon** can not be used directly to solve *infinite-dimensional* constrained optimization problems, e.g., trajectory optimization and optimal control problems. 

## Wrapper



## References
[1] https://sourceforge.net/projects/cocotools/

[2] Dankowicz, H., & Schilder, F. (2013). Recipes for continuation. Society for Industrial and Applied Mathematics.

[3] Li, M., & Dankowicz, H. (2018). Staged construction of adjoints for constrained optimization of integro-differential boundary-value problems. SIAM Journal on Applied Dynamical Systems, 17(2), 1117-1151.

[4] Li, M., & Dankowicz, H. (2020). Optimization with equality and inequality constraints using parameter continuation. Applied Mathematics and Computation, 375, 125058.

[5] https://www.mathworks.com/help/optim/ug/fmincon.html
