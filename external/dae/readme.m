this toolbox combines the ddaecoll_v3 toolbox and the alg_dae_v3 toolbox in DDAE constructors

we modify alg_dae_v3 toolbox such that it generates either inequality or zero (default) function

note that the alg_dae_v3 assumes the output of g function is of the same dimension as y. I used 
a separate id such that we can append several g functions to a ddaecoll segment. You may see the
moon lander example for more details

