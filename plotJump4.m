%plot jumps in scatter plot - jumpSP
clear

load('jumpSP'); 

allRho=[];
filteredRho=[];
allJump=[];
filteredJump=[];

%jump of every merge in index 1. rho in index 2. make one vector with all
%jumps and one without the smallest in each group. 

for i = 2:length(jumpSP)
    jump = jumpSP{i,1};
    rho = jumpSP{i,2};
    
    m = min(jump);
    ind = find(jump==m);
    
    for j = 1:length(jump)
        if rho ~= 49
            allRho = [allRho rho];
            allJump = [allJump jump(j)];
            if j ~= ind && i ~=119
                filteredRho = [filteredRho rho];
                filteredJump = [filteredJump jump(j)];
            end
        end
    end
end

%format plotting
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',24)

hold on
scatter(allRho,allJump, 'r')
scatter(filteredRho, filteredJump, 'b', 'filled')
hold off
xlabel('\rho')
ylabel('rmsd')
title('Pre/Post Merge RMSD')