%Estimate the number of clusters as a function of rho and N by looping over
%rho and doing an exponential fit
clear, clf

load('NcStickyMED12.mat'); %load the data
Nmax=size(NcStickyMED12,1); %Number of data sets to use
x=6:5+Nmax; %N domain
C=zeros(1,50); lambda=zeros(1,50); %storage for Nc parameters

for i=1:50
    y=NcStickyMED12(:,i); %Get data at fixed rho as fn of N
    yfit=fit(x',y,'exp1'); %Fit a 1 term exponential
    C(i)=yfit.a; lambda(i)=yfit.b;  %Get vaslue of C and lambda from fit
end

%Plot C
figure(1)
rho0=4; %rho value to start using data at
rhoEnd=49; %rho value to end at
xfit=linspace(rho0,80,1000);
form='b*x^(-c)'; %power law fit
%form='b*exp(-c*x)';
Cfit=fit((rho0:rhoEnd)',C(rho0:rhoEnd)'-C(rhoEnd-1),form, 'startpoint',[0,0],'lower',[0,0])

loglog(rho0:rhoEnd, C(rho0:rhoEnd),'o', 'linewidth',2)
hold on
loglog(xfit, Cfit(xfit)+C(rhoEnd-1),'linewidth',2)
hold off
xlabel('\rho'); ylabel('C(\rho)'); 
title('Number of Clusters Fit Parameter, C(\rho)');
legend('Empirical Result', 'Fitted Curve','location','east')

%Plot lambda. Try to fit. for lambda>8
figure(2) 
%rho0=3; %rho value to start using data at
xfit=linspace(2,80,1000);
form='a-b*x^(-c)'; %power law fit
Lfit=fit((rho0:rhoEnd)',lambda(rho0:rhoEnd)',form, 'startpoint',[2,1,1.5],'lower',[0,0,0])
hold on
plot(1:50, lambda(1:50), 'o','linewidth',2)
plot(xfit, Lfit(xfit),'linewidth',2)
axis([0,50,0,2])
hold off
xlabel('\rho'); ylabel('\lambda(\rho)');
title('Number of Clusters Fit Parameter, \lambda(\rho)');
legend('Empirical Result', 'Fitted Curve','location','east')

%Plot surface for Nc for rho>=4 using fitted curves
% figure(3)
% N=linspace(6,12,1000);
% r=linspace(4,60,1000);
% Nc=@(x,y) (Cfit(y)+C(49))*exp(Lfit(y)*x);
% fsurf(Nc, [6,12,4,60])