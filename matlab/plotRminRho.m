% plot merge rmin as fn of rho
clear; clf; 

% rho = [4,6,7,8,9,11,20,27,32,35,38,44]; 
% SHS = [1.633,1.4142,1.4142,1.4142,1.2892,1.2737, 1.1785,  1.1269, 1.0886, 1.0866, 1.0853, 1.0833];
% r = [1.6176, 1.3614, 1.342, 1.2563, 1.239, 1.207, 1.143, 1.104, 1.0792, 1.0752, 1.0663, 1.063]; 
% s = [1.619,1.3776,1.36,1.2735,1.24,1.21,1.1405,1.101,1.0782,1.0745,1.0657,1.063];
% w= [1,1,1,1,1,1,10,1,1,10,10,10]; 
% w2= [5,1,1,0,5,0,5,1,1,1,1,10];

rho = [3.5, 5.1, 6.75, 10.4, 24.4, 27.1, 35.26, 38.13];
r = [ 1.58, 1.4, 1.3, 1.2, 1.1, 1.09, 1.075, 1.0645];
m = [3.4, 5.09, 6.75, 9.5, 24.1, 27.5, 35.16, 38.23];
s = [1.57, 1.397, 1.285, 1.203, 1.094, 1.09, 1.074, 1.0639];


f = fit(rho',r','a*x^b+c', 'start', [2,-1,1])% 'weight',w)
%f2 = fit(rho',SHS','a*x^b+c', 'start', [2,-1,1], 'weight',w2)
rp = linspace(3,50,1000); 

% Make figures look better
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',24)

hold on
plot(rho,r,'o')
plot(m,s,'x')
plot(rp, f(rp));
%plot(rp, 2.3*rp.^(-1.05)+1.01)
%plot(rp, 2.3*rp.^(-1.0)+1.01)
%plot(rho, SHS, '+')
%plot(rp,f2(rp));
legend('Morse Measurements', 'Lennard-Jones Measurements', ' 2.1\rho^{-1.05}+1.03') %'SHS r_{min}', ' 2.3\rho^{-0.95}+1.02')
%title('r_{min} Values for a Merge at \rho')
xlabel('\rho'); ylabel('r_{min}'); 
axis([0,50,1,1.8])
hold off
