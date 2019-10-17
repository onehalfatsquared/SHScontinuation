% plot merge rmin as fn of rho
clear; clf; 

rho = [4,6,7,8,9,11,20,27,32,35,38,44]; 
SHS = [1.633,1.4142,1.4142,1.4142,1.2892,1.2737, 1.1785,  1.1269, 1.0886, 1.0866, 1.0853, 1.0833];
r = [1.6176, 1.3614, 1.342, 1.2563, 1.239, 1.207, 1.143, 1.104, 1.0792, 1.0752, 1.0663, 1.063]; 
s = [1.619,1.3776,1.36,1.2735,1.24,1.21,1.1405,1.101,1.0782,1.0745,1.0657,1.063];
w= [1,1,1,1,1,1,10,1,1,10,10,10]; 
w2= [5,1,1,0,5,0,5,1,1,1,1,10];

f = fit(rho',r','a*x^b+c', 'start', [2,-1,1], 'weight',w)
f2 = fit(rho',SHS','a*x^b+c', 'start', [2,-1,1], 'weight',w2)
rp = linspace(3,50,1000); 

% Make figures look better
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',24)

hold on
plot(rho,r,'o')
plot(rho,s,'x')
plot(rp, f(rp));
plot(rho, SHS, '+')
plot(rp,f2(rp));
legend('Morse Measurements', 'Lennard-Jones Measurements', ' 2.3\rho^{-1.05}+1.01', 'SHS r_{min}', ' 2.3\rho^{-0.95}+1.02')
%title('r_{min} Values for a Merge at \rho')
xlabel('\rho'); ylabel('r_{min}'); 
hold off