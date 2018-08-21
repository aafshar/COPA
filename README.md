# COPA: constrained parafac2 for sparse & large datasets
The COPA repository provides the code from our CIKM 2018 paper: "copa: constrained parafac2 for sparse & large datasets", 
By Ardavan Afshar, Ioakeim Perros, Evangelos E. Papalexakis, Elizabeth Searles, Joyce Ho, Jimeng Sun. 

This repository is designed for imposing different constraints on different modes of an irregular tensor (PARAFAC2)

Before running the codes you need to import Tensor Toolbox Version 2.6 which can be downloaded from: https://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html

To start with you need to run: "Run_This.m" file. 
You can apply the following constraints in COPA framework:

1-non-negativity on H, S_k, and V

2-Smoothness on U_k

3-sparsity constraint (l1 and l0) on H, S_k, V

Here is the lists of COPA functions:

create_parafac2_problem: -create a synthetic parafac 2 tensor.

calculate_fit:           -compute the fit for PARAFAC2 tensor input.

claculate_norm:          -compute the norm of a PARAFAC2 tensor.

MSplineBasis:            -This function produce the spline function for subject X_k.

Smooth_COPA:             -Smooth PARAFAC2 where smoothness apply to U_k factor matrix.

fastADMM:                -compute the admm for each mode of a tensor.

COPA                     -This function contains the core part of COPA framework.

COPA_optimizer           -This function is designed to  optimize H,W (S_k) and V.

If you find any bug or error in the codes please send an email to: aafshar8@gatech.edu
