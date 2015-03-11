clear all   % remove all existing variables from memory
clear all   % close all open figures
clc         % clear the MATLAB command window
a = 1;      % define the advection speed
xL = -1;    % define the left interval boundary
xR = +1;    % define the right interval boundary
IMAX = 100; % choose the number of points to discretize the interval width
dx = (xR-xL)/(IMAX -1); % define the mesh spacing
dt = 0.001;  %choose a time step
% Initial contition is given by two piecewise constant
% states, separated by a discontinuity (so-called Riemann problem)
qL = 1;     %left initial state
qR = 0;     %right initial state

for i=1:IMAX
    x(i) = xL + (i-1)*dx; %generate the "mesh" in space
    if(x(i)<=0)
        q(i) = qL;  %assign the left state
    else
        q(i) = qR;  %assign the right state
    end
end

plot(x,q,'o')
time = 0;   %initial time
tend = 0.5; %final time
% time loop (main part of the code)
NMAX = 1000;
for n=1:NMAX
   if ( time+dt > tend)
       dt = tend - time; % reduce the last time step in order to reach tend e???
   end
   if (time>=tend)
       break    % if we have reached tend, the stop the loop
   end
   % evolve the solution to the new time level
   for i=1:IMAX
        if(i==1)
            qnew(i) = qL; % left boundary condition
        elseif(i==IMAX)
            qnew(i) = qR; % right boundary condition
        else
            % explicit central FD scheme
            qnew(i) = q(i) -a*dt/dx*(q(i+1)-q(i-1));
        end
   end
   time = time +dt; % advance time
   q = qnew;        % overwrite the old solution with the new q
   plot(x,q,'o')    % plot the result
   title(sprintf(' Current time = %f ', time))
   drawnow          % force MATLAB to plot the graph immediatly
end
    