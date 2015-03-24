function Q=Energy(T)
global kL kR hc cL cR rhoL rhoR lambdaL lambdaR Tc TL TR
if(T<=Tc)
    % ice
    Q = rhoL*cL*T;
else
    % liquid water
    Q = rhoR*cR*T + rhoR*hc;
end
