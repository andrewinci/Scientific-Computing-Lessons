% Gauss elimination applied to a tridiagonal matrix with the Thomas
% algorithm
% Input:
% a = vector of lower diagonal elements;
% b = vector of the diagonal elements;
% c = vector of upper diagonal elements;
% d = right hand side vector
% Output:
% x = solution of the system

function x=Thomas(a,b,c,d)
N=length(b); % compute the size of the system
% Part I: forward elimination
c(1)=c(1)/b(1);
d(1)=d(1)/b(1);
for i=2:N
   tmp=1/(b(i)-c(i-1)*a(i));
   c(i)=c(i)*tmp;
   d(i)=(d(i)-a(i)*d(i-1))*tmp;
end
% Part II: backward loop
x=zeros(N,1);   % allocate a column vector
x(N)=d(N);
for i=N-1:-1:1
    x(i)=d(i)-c(i)*x(i+1);
end