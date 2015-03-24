function T=Temperature(Q)
global kL kR hc cL cR rhoL rhoR lambdaL lambdaR Tc TL TR

if(Q<0)
    % ice
    T = Q/(rhoL*cL);
elseif(Q>=0 && Q<=rhoR*hc)
    % phase change from ice to water
    T = Tc;
else
    % liquid water
    T = (Q-rhoR*hc)/(rhoR*cR);
end