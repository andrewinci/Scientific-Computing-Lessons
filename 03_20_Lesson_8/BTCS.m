% BTCS scheme for the Stefan problem using the nested Newton technique of
% Casulli and Zanolli

clear all
close all
clc
global kL kR hc cL cR rhoL rhoR lambdaL lambdaR Tc TL TR epsilon
epsilon = 0.1;      % size of the regularization region of the jump

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
gamma=Newton(1);
G1=(Tc-TL)/erf(gamma);
G2=(Tc-TR)/erfc(gamma*sqrt(kL/kR));
s=gamma*2*sqrt(kL*t_end);
xe=linspace(xL,xR,10*IMAX); % mesh for plotting the exact solution
% Set the initial condition
for i=1:IMAX
    T(i) = TR;              % initial temperature of the lake
end

NMAX = 100000;
for n=1:NMAX
    dt = 3600; % 0.45*dx^2/max(kL,kR);
    if(time+dt>t_end)
        dt = t_end-time;     
    end
    if(time>=t_end)
        break
    end
    % compute the temperature from the internal energy
    for i=1:IMAX
        if(T(i)>Tc)
            lambda(i)=lambdaR;
        else
            lambda(i)=lambdaL;
        end
    end
    hold off
    axis([xL xR TL TR])
    plot(x,T,'o')
    title(sprintf('Current time = %f', time))
    drawnow
    
    for i=1:IMAX
        if(i==1)
            lambdap = 0.5*(lambda(i)+lambda(i+1));
            lambdam = lambda(i);
            a(i) = 0;
            b(i) = dt/dx^2*(2*lambdam+lambdap);
            c(i) = -lambdap*dt/dx^2;
            rhs(i) = Q(T(i))+2*lambdam*dt/dx^2*TL;
        elseif(i==IMAX)
            lambdap = lambda(i);
            lambdam = 0.5*(lambda(i)+lambda(i-1));
            a(i) = -lambdam*dt/dx^2;
            b(i) = dt/dx^2*(lambdam+2*lambdap);
            c(i) = 0;
            rhs(i) = Q(T(i))+2*lambdap*dt/dx^2*TR;
        else
            lambdap = 0.5*(lambda(i)+lambda(i+1));
            lambdam = 0.5*(lambda(i)+lambda(i-1));
            a(i) = -lambdam*dt/dx^2;
            b(i) = dt/dx^2*(lambdam+lambdap);
            c(i) = -lambdap*dt/dx^2;
            rhs(i) = Q(T(i));
        end
    end % loop over i
    % We can start with the nested Newton method
    tol = 1e-12*rhoR*hc;    % relative tolerance has to be scaled 
                               % with respect to the number that 
                                  % we are dealing with
    % Initial guess
    T = min(T,Tc-epsilon);
    for iNewton=1:100   % outer Newton iterations
        % The task of outer iterations is only to linearize
        % one of the two non-linear energy functions
        for i=1:IMAX
           % compute the non-linear function f(T) 
           f(i) = Q(T(i))-rhs(i);
           if(i==1)
               f(i) = f(i)+b(i)*T(i)+c(i)*T(i+1);
           elseif(i==IMAX)
               f(i) = f(i)+a(i)*T(i-1)+b(i)*T(i);
           else
               f(i) = f(i)+a(i)*T(i-1)+b(i)*T(i)+c(i)*T(i+1);
           end
        end
        outres = sqrt(sum(f.*f));    % outer residual
        disp(sprintf(' Outer iteration %d, outers = %e', iNewton, outres))
        if(outres<tol)
            break     % if the outer residual is below the tolerance
        end
        Tk = T;  % save the temperature at the outer iteration
        % initial guess for the inner iteration
        T = max(T,Tc-epsilon);
        for inner=1:100
           for i=1:IMAX
              fk(i) = Q1(T(i))-(Q2(Tk(i))+dQ2(Tk(i))*(T(i)-Tk(i)))-rhs(i);
              di(i) = dQ1(T(i))-dQ2(Tk(i));
              if(i==1)
                  fk(i) = fk(i)+b(i)*T(i)+c(i)*T(i+1);
              elseif(i==IMAX)
                  fk(i) = fk(i)+a(i)*T(i-1)+b(i)*T(i);
              else
                  fk(i) = fk(i)+a(i)*T(i-1)+b(i)*T(i)+c(i)*T(i+1);
              end
           end
           inres = sqrt(sum(fk.*fk));
           disp(sprintf('Inner iteration %d, inres = %e', inner, inres))
           if(inres<tol)
               break;
           end
           dT = Thomas(a,b+di,c,fk);    % inner Newton step
           T = T(:)-dT(:);  % update the temperature in Newton loop
        end
    end     % end of the outer iterations    
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
axis([xL xR TL TR])
plot(xe,Te,'r-')    % plot the exact solution
plot(s,Tc,'r+')