clear all   % remove all existing variables from memory
close all   % close all open figures
clc         % clear the MATLAB command windows

a=1;        % define the advection speed
xL=-1;      % define the left interval boundary
xR=+1;      % define the right interval boundary
IMAX=100;   % choose the number of points to discretize the interval with
dx=(xR-xL)/(IMAX-1);    % define the mesh spacing
CFL=0.7;    % choose a CFL number <=1

% Initial condition is given by two piecewise constant states,
% separated by a discontinuity (called Riemann problem)
qL=1;       % left initial state
qR=0;       % right initial state

for i=1:IMAX
    x(i)=xL+(i-1)*dx;   % generate the mesh in space
    if(x(i)<=0)
        q(i)=qL;    % assign the left state
    else
        q(i)=qR;
    end
end

plot(x,q,'o')

time=0;     % initial time
t_end=0.5;  % final time

%time loop (main part of the code)
NMAX=1000;  % choose the number of points to discretize the time interval with
for n=1:NMAX
    dt=CFL*dx/abs(a);   % compute the time step according to the  stability condition
    if(time+dt>t_end)
        dt=t_end-time;  % reduce time step in order to reach t_end exactly
    end
    if(time>=t_end)
        break   % if we have reached t_end, then stop the loop
    end
    % evolve the solution to the new time level
    for i=1:IMAX
        if(i==1)
            q_new(i)=qL;    % left boundary condiion
        elseif(i==IMAX)
            q_new(i)=qR;    % right coundary condition
        else    % Explicit upwind FD scheme
            am=0.5*(a-abs(a));
            ap=0.5*(a+abs(a));
            q_new(i)=q(i)-am*dt/dx*(q(i+1)-q(i))-ap*dt/dx*(q(i)-q(i-1));
        end
    end
    time=time+dt;   % advance time
    q=q_new;    % overwrite the old solution with the new one
    plot(x,q,'o')   % plot the result
    title(sprintf('Current time = %f ', time))
    drawnow     % force MATLAB to plot the graph immediately
end

xe=linspace(xL,xR,10*IMAX); % generate a vector of equidistant points (boundaries included)
for i=1:10*IMAX
    xe(i)=xL+(i-1)*dx/10;   % generate the mesh in space
    if(xe(i)-a*time<=0)     % logical operations
        qe(i)=qL;           % assign the left state
    else
        qe(i)=qR;           % assign the right state
    end
end
hold on     % add another graph in the figure
plot(xe,qe,'r-')    % plot the exact solution
legend('Upwind','Exact')