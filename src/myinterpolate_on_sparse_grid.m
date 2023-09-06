function f_values = myinterpolate_on_sparse_grid(S,Sr,function_on_grid_in,non_grid_points) 
%  interpolate only on inner sparse grids with zero boundary values

%%%% adapted to zero boundary values %%%%
N=size(Sr.knots,1);
nb_pts   = size(non_grid_points,2);
f_values = zeros(nb_pts,1); 

% change to columnwise
function_on_grid_in = function_on_grid_in';
non_grid_points = non_grid_points';

nb_grids=length(S);
global_knot_counter=1;

% loop over the grids
for i=1:nb_grids
    % this is the set of points where I build the tensor lagrange function
    knots_in = S(i).knots_in;
    % I will need the knots in each dimension separately, to collocate the lagrange function.
    % I compute them once for all. As the number of knots is different in each direction, I use a cell array
    knots_per_dim=S(i).knots_per_dim;   % coordinates of knots in each dimension
    % we could just loop on the interpolation points of the current tensor grid, 
    % and for each one evaluate the corresponding lagrangian polynomial on
    % the set of non_grid_points, but this is not convenient. Actually doing this we recompute the same thing over and over! 
    % think e.g. of the lagrangian functions for the knot [x1 y1 z1] and
    % for the knot [x1 y1 z2]: the parts on x and y are the same!
    %
    % Therefore we evaluate each monodim_lagr just once and save all these
    % evaluations in a cell array. We will combine them suitably afterward.
    %
    % Such cell array stores one matrix per direction, i.e. N matrices. 
    % In turn, the n-th matrix contains the evaluations of each monodimensional lagrange polynomial 
    % in direction n on all the n-th coordinates of the non_grid_points.
    % Its dimension will be therefore (number_of_points_to_evaluate) X (number_of_lagrangian_polynomials_in_the_nth_direction)
    %
    % This is actually all the information that we will need to combine to
    % get the final interpolant value.
    %
    mono_lagr_eval=cell(1,N);
    % loop on directions
    for dim=1:N
        K = S(i).m(dim);
        % allocate space for evaluations. Since I will be accessing it one lagrangian polynomial at a time
        % i.e. one knot at a time, it's better to have all information for the same lagrange polynomial
        % on the same column, for speed purposes. Moreover, note that whenever K=1 (one point only in a direction)
        % then the lagrange polynomial is identically one
        if K == 1
            mono_lagr_eval{dim}=ones(nb_pts,1);
        elseif K == 3
            mono_lagr_eval{dim} = 1-non_grid_points(:,dim).^2;
        else
            % loop on each node of the current dimension and evaluate the corresponding monodim lagr polynomial.
            % We will need an auxiliary vector to pick up the current knot (where the lagr pol is centered) and
            % the remaining knots (where the lagr pol is zero). Here we see that mono_lagr_eval it's written
            % one column at a time
            mono_lagr_eval{dim} = zeros(nb_pts,K-2);
            aux=1:K;
            for k=2:K-1  % do not need to compute end-point Lagrange polynomials
                mono_lagr_eval{dim}(:,k-1) = lagr_eval_fast(knots_per_dim{dim}(k),...
                    knots_per_dim{dim}(aux~=k),S(i).m(dim)-1,non_grid_points(:,dim),[nb_pts 1]);
            end
        end
    end
    % now put everything together. We have to take the tensor product of
    % each of the monodim lagr pol we have evaluated. That is, we have to
    % pick one column for each matrix in the cell array and dot-multiply them.
    % all the possible combinations have to be generated !
    %
    % once this is done, we have the evaluation of each multidim lagr
    % polynomial on the non_grid_points, which we will then multiply by the
    % corresponding nodal value and eventually sum everything up.
    %
    % We start by generating the combination of column we need to take. We actually don't
    % need to generate them, but only to recover it from the matrix knots,
    % which already contains all the points of the grid, i.e. all the
    % combinations of 1D points!
    %
    % Given a matrix of points like
    % knots=[a1 a1 b1 b1 a1 a1 b1 b1 ...
    %        a2 b2 a2 b2 .....
    %        
    % combi is
    % combi=[1 1 2 2 1 1 2 2 ...
    %        1 2 1 2 ......
    % again, we exploit the fact that the minimum entry of combi is 1 and that for many directions
    % there is only one point, so if we init combi with ones we're good in many cases
    combi = ones(N,S(i).size_in);
    % the easiest way to recover combi from knots is to proceed one dimension at a time,
    % and mark with a different label (1,2,...K) all the equal points. We need of course as many labels
    % as the number of different points in each dir!
    
    for dim=1:N
        % this is how many points per direction
        K=S(i).m(dim);
        %%% we start from a row of zeroes and we place 1....K in the right
        % positions by summations (each element of the row will be written
        % only once!). Since 1 are already in place, we proceed to place 2 and higher, but only if needed
        if K>1
            for k=2:K-1
                % here we add to the row of "1" either 0 or (k-1) so we get k where needed
                combi(dim,:) = combi(dim,:)+ (k-1)*( knots_in(dim,:)==knots_per_dim{dim}(k) ); 
            end
        end
    end
    combi = combi-1; combi(combi == 0) = 1;
    %%% Now we can do the dot-multiplications among rows, the
    % multiplication by nodal values and the final sum! We proceed one
    % knot at a time
    coeff=S(i).coeff;
    
    for kk=1:S(i).size_in
        f_loc=mono_lagr_eval{1}(:,combi(1,kk));
        for dim=2:N
            f_loc=f_loc.*mono_lagr_eval{dim}(:,combi(dim,kk));
        end
        % recover F, the corresponding value for the interpolating function in function_on_grid, with the global counter
        position = Sr.n_in(global_knot_counter);
        F_value = function_on_grid_in(position,:);

        f_values = f_values + coeff*f_loc*F_value; 
        global_knot_counter=global_knot_counter+1;
    end
end % of for loop on tensor grid
% finally, transpose to comply with output orientation
f_values=f_values';



