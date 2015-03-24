clear all
close all
clc
global kL kR hc cL cR rhoL rhoR lambdaL lambdaR Tc TL TR

% Physical parameters of water (all in SI units)
lambdaL = 2.09;     % heat conductivity of ice
lambdaR = 0.6;      % heat conductivity of water
hc      = 334e3;    % latent heat
rhoL    = 917;      % density of ice
rhoR    = 1000;     % density of water
cL      = 2108;     % heat capacity of ice
cR      = 4187;     % heat capacity of water
kL      = lambdaL/(rhoL*cL);
kR      = lambdaR/(rhoR*cR);
Tc      = 0;        % critial (melting) temperature
TL      = -10;      % air temperature (TL<Tc)
TR      = +5;       % initial lake temperature (TR>Tc)

% Domain
xL      = 0;        % lake surface
xR      = 2;        % bottom of the lake
IMAX    = 100;      % number of control volumes
dx      = (xR-xL)/IMAX; % mesh spacing
x       = linspace(xL+dx/2,xR-dx/2,IMAX);
day     = 24*3600;  % number of seconds in a day

% Plot the exact solution after a certain time
time  = 0;      % initial time
t_end = 2*day;  % final time
gamma=Newton(1)
G1=(Tc-TL)/erf(gamma);
G2=(Tc-TR)/erfc(gamma*sqrt(kL/kR));
s=gamma*2*sqrt(kL*t_end);
xe=linspace(xL,xR,10*IMAX); % mesh for plotting the exact solution
% set the initial condition
for i=1:IMAX
    T(i) = TR;              % initial temperature of the lake
    Q(i) = Energy(T(i));    % initial internal energy of the lake
end

NMAX = 10000;
for n=1:NMAX
    dt = 0.45*dx^2/max(kL,kR);
    if(time+dt>t_end)
        dt = t_end-time;     
    end
    if(time>=t_end)
        break
    end
    % compute the temperature from the internal energy
    for i=1:IMAX
        T(i)= Temperature(Q(i));
        if(T(i)>Tc)
            lambda(i)=lambdaR;
        else
            lambda(i)=lambdaL;
        end
    end
    hold off
    plot(x,T,'o')
    title(sprintf('Current time = %f', time))
    axis([xL xR -10 5])
    drawnow
    for i=1:IMAX
        if(i==1)
            lambdap = 0.5*(lambda(i)+lambda(i+1));
            lambdam = lambda(i);
            fp = lambdap*(T(i+1)-T(i))/dx;
            fm = lambdam*(T(i)-TL)/(dx/2);
        elseif(i==IMAX)
            lambdap = lambda(i);
            lambdam = 0.5*(lambda(i)+lambda(i-1));
            fp = lambdap*(TR-T(i))/(dx/2);
            fm = lambdam*(T(i)-T(i-1))/dx;
        else
            lambdap = 0.5*(lambda(i)+lambda(i+1));
            lambdam = 0.5*(lambda(i)+lambda(i-1));
            fp = lambdap*(T(i+1)-T(i))/dx;
            fm = lambdam*(T(i)-T(i-1))/dx;
        end
        % conservative finite volume update of the internal energy
        Q(i) = Q(i) + dt/dx*(fp-fm);
    end % loop over i
    % advance in time
    time = time + dt;
end % time loop

for i=1:10*IMAX
   if(xe(i)<s)
       Te(i)=TL+G1*erf(xe(i)/(2*sqrt(kL*time)));
   else
       Te(i)=TR+G2*erfc(xe(i)/(2*sqrt(kR*time)));
   end
end
hold on
plot(xe,Te,'r-')    % plot the exact solution
plot(s,Tc,'r+')
