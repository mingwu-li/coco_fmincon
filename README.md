# COCO wrapper functions for fmincon in MATLAB

## COCO
Continuation Core (COCO) is a MATLAB-based open-source package for computational nonlinear analysis of dynamical systems. It uses atlas algorithms to cover the solution manifolds of nonlinear systems. It supports the detection of bifurcation points on solution manifolds and the switch from one solution branch to the other.

A unique feature of COCO is the embedded construction philosphy of nonlinear systems, where a large problem is assembed by small subproblems with weak couplings. This object-oriented construction paradigm enables us to build a composite problem from building blocks.

COCO provides a predefined library for such building blocks. Specifically, it has constructors for collocation problems (coll), boundary value problems (bvp), periodic orbit problems (po) and hybrid system periodic orbit problems (hspo). Further, COCO provides constructors for the adjoint operators of these problems. Such adjoints can be used to formulate the first-order necessary conditions to constrained optimization problems. Please refer [1-3] for more info of COCO.

The necessary conditions are often given in the form of bvps, which can be solved using successive continuation methods [3-4]. Here, we provide an alternative for solving constrained optimization problems. COCO constructors are used to generate (discretized) constrained optimization problems, which are then solved using fmincon in MATLAB.


## fmincon
**fmincon** is a MATLAB function for finding the minimum of *finite-dimensional* constrained optimization problem. It is a nonlinear programming solver. More details about **fmincon** can be found at MATLAB documentation [5]. Note that **fmincon** can not be used directly to solve *infinite-dimensional* constrained optimization problems, e.g., trajectory optimization and optimal control problems.

## Wrapper functions
An optimization problem includes objective function(al), dimensional equality constraints, and inequality constraints. The *zero* functions in COCO can be used to represent equality constraints and *inequality* functions can be used for inequality constraints. As for objective function(al), one can use monitor functions to represent it. 

Based on these observations, we write wrapper functions for fmincon. Given an problem constructed using COCO, all zero/inequality functions are interpreted as equality/inequality constraints, and a monitor function with identifyer is used to define the objective.

A few examples including a bang-bang optimal control problem are included in this repo. More info about the wrapper functions and examples can be found at doc folder.



## References
[1] https://sourceforge.net/projects/cocotools/

[2] Dankowicz, H., & Schilder, F. (2013). Recipes for continuation. Society for Industrial and Applied Mathematics.

[3] Li, M., & Dankowicz, H. (2018). Staged construction of adjoints for constrained optimization of integro-differential boundary-value problems. SIAM Journal on Applied Dynamical Systems, 17(2), 1117-1151.

[4] Li, M., & Dankowicz, H. (2020). Optimization with equality and inequality constraints using parameter continuation. Applied Mathematics and Computation, 375, 125058.

[5] https://www.mathworks.com/help/optim/ug/fmincon.html
