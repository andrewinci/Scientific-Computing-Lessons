% Matrix-free conjugate gradient method for the soluzion of the linear system Ax=b
% Input:
%   A = system matrix
%   b = right hand side
% Output
%   x = solution vector
% The user need to provide the matrix-vector product as a function "matop"

function x=CGop2D(b)
N=numel(b);             % compute the number of unknowns
x=b;                    % initial guess
r=b-matop2D(x);         % initial residual
p=r;                    % initial search direction
alpha=sum(sum(r.*r));   % square of the norm
for k=0:N-1
    if(sqrt(alpha)<1e-14)
        return              % if the norm is small enough stop iterating
    end
    v=matop2D(p);                     % compute the matrix-vector product
    lambda=alpha/sum(sum(p.*v));    % 1D minimization
    x=x+lambda*p;                   % compute the new iterate
    r=r-lambda*v;                   % compute the new residual
    alpha_k=alpha;                  % save the old alpha
    alpha=sum(sum(r.*r));           % compute the new alpha
    p=r+alpha/alpha_k*p;            % compute the new search direction
end