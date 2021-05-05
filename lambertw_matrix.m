% Written by Shiming Duan, July. 31, 2010
% Dept. Mechanical Engineering, University of Michigan  
% Email: duansm@umich.edu                               
% For Matlab Toolbox of Time-Delay System Supplements 
% K is the branch ID and is a scalar. If X is a matrix and has full rank,
% K is the same for each eigenvalue of X. If X is not invertible, for its zero
% eigenvalues, 0 will be used for the branch ID for this branch. 

function W = lambertw_matrix(K,X)

if nargin == 1 % Only calculate principal branches
    X = K;
    K = 0;
end

rank_X= rank(X); 
dim_X = length(X); % system order
LMW = zeros(dim_X);
[V,D] =  jordan(X); % transform to Jordan form

% check superdiagonal terms
super_diag = [];
if dim_X==1;
    super_diag = 0;
else
    for i = 2:dim_X
        super_diag = [super_diag D(i-1,i)];
    end
end

% calculate the Lambert W func. for diagnal terms
for i = 1:dim_X
    if D(i,i)==0 
        % use principal branch for 0 eigenvalue
        LMW(i,i)=lambertw(0,D(i,i)); % warning('Ad is singular, hybrid branch is used')
    else
        LMW(i,i)=lambertw(K,D(i,i)); 
    end
end

% calculate the Lambert W func. for other upper triangle terms
if sum(super_diag)==0
    flag =0; % all jordan blocks have a size of 1
else
    flag =1; % some jordan blocks have a size > 1
    block_start = 1; % where an jordan block begins
    block_end = 1; % where an jordan block ends
    for j=2:dim_X
        if j>=block_end
            if super_diag(j-1)==1
                block_start = j-1;
                if ~isempty(find(super_diag(j-1:end)==0))
                    block_end = min(find(super_diag(j-1:end)==0))+j-2;
                else
                    block_end = dim_X;
                end
                for k1 =block_start:block_end
                    for k2 = k1: block_end
                        if D(block_start,block_start)~=0
                            LMW(k1,k2)=der_lambertw(K,D(block_start,block_start),k2-k1)/factorial(k2-k1); % (k2-k1)th deritive of lambertw(D(i))/(k2-k1)! 
                        end
                    end
                end
            end
        end
    end
end

W = V*LMW*V^-1; % similarity transformation



  