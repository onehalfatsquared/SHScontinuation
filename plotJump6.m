%plot jumps in scatter plot - jumpSP
clear

load('jumpSP6'); 

allRho=[];
filteredRho=[];
allJump=[];
filteredJump=[];

%get rid of coarse values
% for i = 148-42:148
%     jumpSP6{i,:}=;
% end
    

%jump of every merge in index 1. rho in index 2. make one vector with all
%jumps and one without the smallest in each group. 

for i = 1:length(jumpSP6)
    jump = jumpSP6{i,1};
    rho = jumpSP6{i,2};
    
    m = min(jump);
    ind = find(jump==m);
    
    for j = 1:length(jump) 
        if rho ~= 49 && all(i~=148-48:148) && i~=190 && i~=30 && all(i~=192:194) && all(i~=197:199) && i~=201
            allRho = [allRho rho];
            allJump = [allJump jump(j)];
            if j ~= ind && i~=118 && all(i~=148-42:148)
                filteredRho = [filteredRho rho];
                filteredJump = [filteredJump jump(j)];
                i;
                if jump(j)<0.01
                    i
                end
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
%title('Pre/Post Merge RMSD')

%plot log scale
figure(2)
hold on
scatter(log(allRho),log(allJump), 'r')
scatter(log(filteredRho), log(filteredJump), 'b', 'filled')
hold off
xlabel('log(\rho)')
ylabel('log(rmsd)')