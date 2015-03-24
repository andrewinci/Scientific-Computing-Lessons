% internal energy as a function of temperature using a piecewise linear
% regularization in a small interval [Tc-epsilon, Tc+epsilon]

function y=Q(T)
global kL kR hc cL cR rhoL rhoR lambdaL lambdaR Tc TL TR epsilon

if(T<Tc-epsilon)
    % ice
    y = rhoL*cL*(T-Tc);
elseif(T>=Tc+epsilon)
    % water
    y = rhoR*cR*(T-Tc)+rhoR*hc;
else
    % phase transition in the interval [Tc-epsilon, Tc+epsilon]
    dQdT = ((rhoR*cR*(Tc+epsilon-Tc)+rhoR*hc)-(rhoL*cL*(Tc-epsilon-Tc)))/(2*epsilon);
    y = (rhoL*cL*(Tc-epsilon-Tc))+dQdT*(T-(Tc-epsilon));
end
