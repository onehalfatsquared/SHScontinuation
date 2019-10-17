function [kappa] = stickyEval(E,rho)
%Evaluate the sticky parameter for Morse Potential as function of rho and E
kappa=exp(E)/(sqrt(2*rho^2*E));
end

