% Newton method with a simple globalization strategy 
% for the solution of the equation g(gamma)=0

function gamma=Newton(gamma0)
gamma=gamma0;   % initial guess
tol=1e-10;      % tolerance
for iNewton=1:100
   gk=g(gamma); % function at the current Newton iterate
   res=abs(gk); % compute the residual
   disp(sprintf(' iter = %d, current residual = %e', iNewton, res))
   if(res<tol)
      break     % we have found a solution 
   end
   dgamma=-gk/dg(gamma);    % compute the Newton step
   delta=1;                 % first, try an entire Newton step
   while(abs(g(gamma+delta*dgamma))>res)
      delta=delta/2;        % half the step size while the residual increases
   end
   gamma=gamma+delta*dgamma;    % now do a modified Newton step
end