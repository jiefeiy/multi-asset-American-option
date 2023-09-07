# multi-asset-American-option

This repo includes the algorithms, figures, tables of the paper 'On Sparse Grid Interpolation for American Option Pricing with
Multiple Underlying Assets' by Jiefei Yang and Guanglian Li. 

The main goal is to price American arithmetic and geometric put options with d assets. 

The folder 'algorithms' contains implementation of Algorithm 1 of the paper with various quadrature methods including 

1.  Genz-Keister sparse grid quadrature,
2.  Gauss-Hermite sparse grid quadrature,
3.  normal Leja sparse grid quadrature,
4.  QMC with scramble Sobol sequence,
5.  preintegration strategy and QMC with scramble Sobol sequence.

Examples of using these algorithms are included in the folder 'examples'.

The following two Matlab Toolboxes are needed. 

1. "Parrallel Computing Toolbox",
2. "Statistics and Machine Learning Toolbox".

References:
The implementation of sparse grids is based on the Sparse Grids Matlab Kit, 
a MATLAB toolbox for high-dimensional quadrature and interpolation, see https://sites.google.com/view/sparse-grids-kit.
