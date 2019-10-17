%Scatter plot of rmsd values before and after a merge as fn of rho.
clear

%data
rho = [4,4,6,6,6,35,35,38,38,38,38,38,38];
rmsd = [0.2552, 0.2685, 0.1117, 0.112, 0.109, 0.0218, 0.0228, 0.0202, 0.0194, 0.0195, 0.0192, 0.0201, 0.0194];

%format plotting
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',24)

%plot
scatter(rho,rmsd)
xlabel('\rho')
ylabel('rmsd')