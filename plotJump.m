%plot jumps in scatter plot
clear

load('jump'); 

x=[];
y=[];

c = [ones(1,25)*0.04 ones(1,25)*0.01];

for i=1:50
    J = jump{i, 2};
    for j = 1:length(J)
        if J(j)>c(i)
            x=[x, i];
            y=[y,J(j)];
        end
    end
end

%format plotting
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',24)

scatter(x,y)
xlabel('\rho')
ylabel('rmsd')
title('Pre/Post Merge RMSD')