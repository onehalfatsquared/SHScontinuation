%plot morse potential for various values of rho and E
clear

r = linspace(0.95,1.05,1000);
rho = [10,20,30,40,50];
kappa = 40; 
E = zeros(1,length(rho));
E = rho*2;
% for i = 1:length(rho)
%     E(i)=stickyNewton(8, rho(i), kappa);
% end

hold on
for i =1:length(rho)
    M = E(i)*(exp(-2*rho(i)*(r-1))-2*exp(-rho(i)*(r-1)));
    plot(r,M)
end

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultAxesFontSize',20)

axis([0.95,1.05,-100,100])
legend('\rho = 10','\rho = 20','\rho = 30','\rho = 40','\rho = 50');
xlabel('r'); ylabel('U(r)');
title('Morse Potential')