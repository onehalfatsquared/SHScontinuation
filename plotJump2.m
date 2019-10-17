%plot jumps in scatter plot - jump2
clear

load('jump2'); 

x=[];
y=[];

c = zeros(1,50); c(49)=10;
%jump2{49,3}=jump2{49,3}*10;
tol=1e-3;

for i=1:50
    J = jump2{i, 2};
    r = jump2{i,3}; 
    if isempty(r)
        c(i)=0;
    else
        if max(r)-min(r)>tol
            c(i) = min(r) + 1e-2;
        else
            c(i) = min(r) -1e-2;
        end
    end
    c(49)=10;
    for j = 1:length(J)
        if r(j)>c(i)
            x=[x, i];
            y=[y,J(j)];
            if J(j)<0.01
                i,j
            end
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