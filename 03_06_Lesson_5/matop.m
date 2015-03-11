% Matrix-vector product operator for the matrix-free version of the
% CG-method
% Input: T = temperature
% Output: AT = matrix-vector product = left hand side of the scheme

function AT=matop(T)
global dt dx b IMAX     % define some useful global variables

AT=T;
for i=1:IMAX
   if(i==1)                 % left boundary
       bp=0.5*(b(i)+b(i+1));    % average b on the right
       fp=bp*(T(i+1)-T(i))/dx;  % right flux
       bm=b(i);                 % b on the left
       fm=bm*(T(i)-0)/(dx/2);   % left flux
   elseif(i==IMAX)          % right boundary
       bp=b(i);                 % b on the right
       fp=bp*(0-T(i))/(dx/2);   % right flux
       bm=0.5*(b(i)+b(i-1));    % average b on the left
       fm=bm*(T(i)-T(i-1))/dx;  % left flux
   else
       bp=0.5*(b(i)+b(i+1));    % average b on the right
       fp=bp*(T(i+1)-T(i))/dx;  % right flux
       bm=0.5*(b(i)+b(i-1));    % average b on the left
       fm=bm*(T(i)-T(i-1))/dx;  % left flux
   end
   AT(i)=AT(i)-dt/dx*(fp-fm);
end