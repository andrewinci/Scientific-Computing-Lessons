% MATLAB function for Gauss elimination with pivoting
% Input:
% A =  system matrix
% b = right hand side vector
% Output:
% x = solution of the system
% Function name: Gauss
% File name: Gauss.m

function x=Gauss(A,b)
N=length(b);    % get the number of unknowns
C=[A,b];        % [] is an array constructor

% Part 1: Forward elimination (produces an upper triangular matrix)
for i=1:N
    maxrow=0;
    maxval=0;
    for j:1=N
       if(abs(C(j,i))>maxval)
           maxval=abs(C(j,i));
           maxrow=j;
       end
    end
    if(maxval==0)
        disp('Matrix is singular')
        x=[];   % the empty set
        return
    else        % exchange of the rows i,j
        j=maxrow;
        temp=C(j,:);
        C(j,:)=C(i,:);
        C(i,:)=temp;
    end
    % divide row i by the diagonal element
    C(i,:)=C(i,:)/C(i,i);
    % loop over all rows j>i
    for j=i+1:N
        C(j,:)=C(j,:)-C(j,i)*C(i,:);    % eliminate all elements in row i
    end
end
% Now the matrix C is in upper triangular form
% Part 2: back substitution
for i=N:-1:1    % start in row N and loop backward to 1
    for j=i-1:-1:1
       C(j,:) =C(j,:)-C(j,i)*C(i,:);
    end
end
% The solution is now contained in the right part of C
x=C(:,N+1);