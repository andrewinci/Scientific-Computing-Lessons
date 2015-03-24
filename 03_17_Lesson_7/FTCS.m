% FTCS method for the linea heat conduction equation
clear all           % removes all variables from memory
close all           % closes all figures
clc                 % clears the command windows
xL=-1;              % left boundary
xR=1;               % right boundary
IMAX=100;           % number of spatial control volumes
dx=(xR-xL)/IMAX;    % mesh spacing
x=linspace(xL+dx/2, xR-dx/2, IMAX);  % generate the grid
%b=1; % b = lambda/(rho*c_v)
bL=10; % left coefficient 
bR=1; % right coefficient
TL=100; % left temperature
TR=20; % right temperature
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
plot(x,T,'o') % plot the initial condition

time=0; % initial time
tend=70.0; % final time
NMAX=100000; % max number of time steps
for i=1:NMAX
    % choose a stable time step
    dt=0.45*dx^2/max(b);
    if(time+dt>tend)
        dt=tend-time; % reduce dt to reach tend exactly
    end
    if(time>=tend)
        break
    end
    for i=1:IMAX
        if(i==1)
            % left boundary
            bp=0.5*(b(i)+b(i+1));
            fp=-bp*(T(i+1)-T(i))/dx;
            %fm=0; % zero heat flux (Neumann BC)
            bm=b(i);
            fm=-bm*(T(i)-TL)/(dx/2); % impose T (Dirichlet BC)
        elseif(i==IMAX)
            % right boundary
            %fp=0; % zero heat flux (Neumann BC)
            bp=b(i);
            fp=-bp*(TR-T(i))/(dx/2); % impose T (Dirichlet BC)
            bm=0.5*(b(i)+b(i-1));
            fm=-bm*(T(i)-T(i-1))/dx;
        else
            % inside the domain
            bp=0.5*(b(i+1)+b(i));
            fp=-bp*(T(i+1)-T(i))/dx;
            bm=0.5*(b(i)+b(i-1));
            fm=-bm*(T(i)-T(i-1))/dx;
        end
        % Finite volume scheme
        Tnew(i)=T(i)-dt/dx*(fp-fm);
    end
    % overwrite the old solution with the new one
    T=Tnew;
    % advance time
    time=time+dt;
    plot(x,T,'o')
    title(sprintf('Current time = %f', time))
    drawnow
end