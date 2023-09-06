function pattern = generate_pattern_in(m,sz_in)

% pattern = generate_pattern(m)
%
% given m=[m1 m2 m3 m4 ... mN] generares a matrix that can be used to generate the cartesian product
% of {1,2,...,m1} x {1,2,...,m2} x {1,2,...m3} x ....
% 
% e.g.
%
% generate_pattern([3 2 2])
% 
% pattern =
% 
%       1      2      3      1      2      3      1      2      3      1      2      3
%       1      1      1      2      2      2      1      1      1      2      2      2
%       1      1      1      1      1      1      2      2      2      2      2      2


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2022 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



N=length(m);

% it is convenient from a computational point of view to generate the pattern as unsigned int to create the pattern
if max(m)<=intmax('uint8')
    pattern=zeros(N,sz_in,'uint8');
elseif max(m) <=intmax('uint16')
    warning('SparseGKit:uint16','more than 255 points are asked in one direction, using uint16 to handle this')
    pattern=zeros(N,sz_in,'uint16');
else
    warning('SparseGKit:double','more than 65535 points are asked in one direction, using double precision to handle this')
    pattern=zeros(N,sz_in);
end


% the algorithm works one direction at a time. at the n-th iteration the n-th row of pattern is generated,
% by repeating q times the vector BASE, which containes itselt a sequence,
% obtained repeating p times each number from j=1 to j=m(n), e.g. 
% generate_pattern([3 2 2])
% 
% pattern =
% 
%       1      2      3      1      2      3      1      2      3      1      2      3
%       1      1      1      2      2      2      1      1      1      2      2      2
%       1      1      1      1      1      1      2      2      2      2      2      2
mm = m-2; 
mm(mm<0) = 1;

for dim=1:N
    p = prod([1 mm(1:dim-1)]);
    q = prod([mm(dim+1:end) 1]);
    
    % length of base vector
    lb=p*mm(dim);
    base = zeros(1,lb);
    
    % generate base vector
    bb=1;
    if m(dim)<2
        base(bb:bb+p-1) = 1;
    else 
        for j=2:m(dim)-1
            base(bb:bb+p-1)=j;
            bb=bb+p;
        end
    end
    
    % repeat base vector
    pp=1;
    for j=1:q
        pattern(dim,pp:pp+lb-1)=base;
        pp=pp+lb;
    end
    
end