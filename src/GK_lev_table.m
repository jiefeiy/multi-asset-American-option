function A = GK_lev_table(rows,cols)

% A = GK_lev_table(rows,cols)
%
% In the tabulated GK knots, many levels have the same knots (with weights identical 
% up to 10th-11th digit. For each group of levels with same number of knots we
% select 1 representative, and define a map as follows:
% 
% MAP=
% i:   l:     nb_knots:   order: 
% 0     0     0             0
% 1     1     1             1
% 2     3     3             3
% 3     8     9             15 
% 4    15    19             29
% 5    25    35             51
%
% this function returns the (row,col) submatrix of M


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2022 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


M=[
0     0     0   0
1     1     1   1
2     3     3   3
3     8     9   15
4    15    19   29
5    25    35   51 ];

A=M(rows,cols);

end
