# multi-asset-American-option

This repository includes the algorithms, figures, tables of the paper 'On Sparse Grid Interpolation for American Option Pricing with Multiple Underlying Assets' by Jiefei Yang and Guanglian Li. 
The code is published to ensure the reproducibility of the numerical examples. 

 - The algorithm can be used to price American arithmetic and geometric put options with d assets.
 - The following two Matlab Toolboxes are needed. 
    * "Parrallel Computing Toolbox",
    * "Statistics and Machine Learning Toolbox". 
 - Before running examples, you will need to add necessary files to the MATLAB path. Just execute
```
init.m
```
## Main algorithms 

The folder `./algorithms` contains implementation of Algorithm 1 of the paper with various quadrature methods including 

1.  Genz-Keister sparse grid quadrature,
2.  Gauss-Hermite sparse grid quadrature,
3.  normal Leja sparse grid quadrature,
4.  QMC with scramble Sobol sequence,
5.  preintegration strategy and QMC with scramble Sobol sequence.

## Quick start

To see a 2-d example with the above algorithms, run `./examples/examples_diff_quad.m`. The following is expected to display in the command window:
```
The results of pricing values are 
1. The value using Genz-Keister sparse grid quadrature is 3.188
2. The value using Gauss-Hermite sparse grid quadrature is 3.1825
3. The value using normal Leja sparse grid quadrature is 3.1833
4. The value using QMC with scramble Sobol sequence is 3.1818
5. The value using preintegration strategy and QMC with scramble Sobol sequence is 3.1877
```
(The results of 4,5 might be slightly different due to the randomness.)

## Reproducibility of the numerical examples

To get the results in Section 5 (Numerical experiments) of the paper, run executable scripts in `./numerical_tests`.

- For Section 5.1, get into the directory `./numerical_tests/section5_1/` first, then 
    * run `t1_arithm_bask_put.m` to obtain values in Table 1 of the paper;
    * run `t2_geo_bask_put.m` to obtain values in Table 2 of the paper;
    * run `t3_exact_geo_bask_put.m` to obtain the reference Bermudan price in Table 2 of the paper;
    * run `tab1_results_arithBask.m` and `tab2_results_geoBask.m` to compute relative errors of Table 1 and 2, respectively.

- For Section 5.2, get into the directory `./numerical_tests/section5_2/` first, then 
    * run `convergence_plot.m` to obtain Figure 4 of the paper.

- For Section 5.3, get into the directory `./numerical_tests/section5_3/` first, then 
    * run `t1_fig_5a_onestep_error.m` to obtain Figure 5(a);
    * run `t2_fig_5b_total_error.m` to obtain Figure 5(b);
    * run `t3_fig_6_convergence_rqmc.m` to obtain Figure 6. 

- For Section 5.4, get into the directory `./numerical_tests/section5_4/` first, then 
    * run `t1_different_beta.m` to obtain Figure 7(a);
    * run `t2_different_L.m` to obtain Figure 7(b);
    * run `t3_different_K.m` to obtain Figure 7(c).

## Reproducibility of Figure 1,2,3
First, get into the directory `./figures123/`, then 
 * run `fig1_b_and_c.m` to obtain Figure 1(b)(c).
 * run `fig2_transform_gridpoints.m` to obtain Figure 2(a)(b).
 * run `fig3_compare_num_points.m` to obtain Figure 3.
 * To obtain Figure 1(a), the toolbox "Algorithm 847: Spinterp: piecewise multilinear hierarchical sparse grid interpolation in MATLAB" (https://dl.acm.org/doi/abs/10.1145/1114268.1114275) is required. After installation, run `fig1_a_maxnorm_no_boundary.m`.

## References
The implementation of sparse grids is based on the Sparse Grids Matlab Kit, 
a MATLAB toolbox for high-dimensional quadrature and interpolation, see https://sites.google.com/view/sparse-grids-kit.
