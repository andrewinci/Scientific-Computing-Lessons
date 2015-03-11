% BTCS method for the linea heat conduction equation
clear all           % removes all variables from memory
close all           % closes all figures
clc                 % clears the command windows
xL=-1;              % left boundary
xR=1;               % right boundary
IMAX=100;           % number of spatial control volumes
dx=(xR-xL)/IMAX;    % mesh spacing
x=linspace(xL+dx/2, xR-dx/2, IMAX);  % generate the grid
%b=1;               % b = lambda/(rho*c_v)
bL=10;              % left coefficient 
bR=1;               % right coefficient
TL=100;             % left temperature
TR=20;              % right temperature
% initial condition
for i=1:IMAX
   if(x(i)<=0) 
       T(i)=TL;
       b(i)=bL;
   else
       T(i)=TR;
       b(i)=bR;
   end
end
plot(x,T,'o')       % plot the initial condition

time=0;             % initial time
tend=7000.0;           % final time
NMAX=100000;        % max number of time steps
for i=1:NMAX
   % choose a stable time step
   dt=0.45*dx^2/max(b);
   if(time+dt>tend)
       dt=tend-time;    % reduce dt to reach tend exactly
   end
   if(time>=tend)
       break
   end
   
   % assemble the coefficients of the matrix A of the linear system
   
   for i=1:IMAX
      if(i==1)
          % left boundary
          bp=0.5*(b(i)+b(i+1));
          bm=b(i);
          di(i)=1+(2*bm+bp)*dt/dx^2;        % diagonal
          up(i)=-bp*dt/dx^2;                % upper diagonal
          rhs(i)=T(i)+2*bm*dt/dx^2*TL;      % boundary condition
      elseif(i==IMAX)
          % right boundary
          bp=b(i);
          bm=0.5*(b(i)+b(i-1));
          lo(i)=-bm*dt/dx^2;                % lower diagonal
          di(i)=1+(bm+2*bp)*dt/dx^2;        % diagonal
          up(i)=0;
          rhs(i)=T(i)+2*bp*dt/dx^2*TR;      % right hand side
      else
          % inside the domain
          bp=0.5*(b(i+1)+b(i));
          bm=0.5*(b(i)+b(i-1));
          lo(i)=-bm*dt/dx^2;                % lower diagonal
          di(i)=1+(bm+bp)*dt/dx^2;          % diagonal
          up(i)=-bp*dt/dx^2;                % upper diagonal
          rhs(i)=T(i);                      % right hand side
      end
   end
   % Finite volume scheme
   T=Thomas(lo,di,up,rhs);                  % solve the linear system
   % advance time
   time=time+dt;
   plot(x,T,'o')
   title(sprintf('Current time = %f', time))
   drawnow
end