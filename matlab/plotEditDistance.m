%plot edit distance between two graphs as a function of range parameter
clear

%Set up total num of clusters for normalization
num_clusters(6)=2; num_clusters(7)=5; num_clusters(8)=13; num_clusters(9)=52;
num_clusters(10)=263; num_clusters(11)=1659;

%set up which edit distance we are computing
type = "R"; %Either MM ML or LL or R for comparisons. Morse and Lennard-Jones, random for comp
n = 8; %Number of particles
N = num_clusters(n); %normalizing factor
types = ["MM", "LL", "ML"];
count = 1;
d = zeros(9,50); 

for test = 1:3
    type = types(test);
    d1 = comparePartitions(n, 1, 2, type); d1 = d1/N; %low med
    d(count,:) = d1; count = count+1;
    d2 = comparePartitions(n, 1, 3, type); d2 = d2/N; %low high
    d(count,:) = d2; count = count+1;
    d3 = comparePartitions(n, 2, 3, type); d3 = d3/N; %med high
    d(count,:) = d3; count = count+1;
end

d4 = comparePartitions(n, 2, 2, "R"); d4 = d4/N;

%format plotting
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',24)

%normalize and plot as fn of rho
figure(1)
hold on
M = max(max(d));
count = 1;
ccount = 0; 
x = 1./(1:50);
col = ["blue", "g", "r"];
ls = ["-", "--", ":"];

for test1 = 1:3
    for test = 1:3
        ccount = rem(count,3)+1;
        plot(x,d(count,:), 'color', col(test1), 'linestyle', ls(ccount));
        count = count+1;
    end
end

plot(x,d4, 'k')
M = max(max(M,d4));

xlabel('\rho^{-1}')
ylabel('Edit Distance')
title('N=10')
%legend('Low to Medium', 'Low to High', 'Medium to High');
axis([0 1 0 M+1])
xticks(0:0.1:1);
xticklabels({'0','','0.2','','0.4','','0.6','','0.8','','1'})
hold off

%normalize and plot as fn of log(rho)
figure(2)
hold on
M = max(max(d));
count = 1;
ccount = 0; 
x = log(1:50);
col = ["blue", "g", "r"];
ls = ["-", "--", ":"];

for test1 = 1:3
    for test = 1:3
        ccount = rem(count,3)+1;
        plot(x,d(count,:), 'color', col(test1), 'linestyle', ls(ccount));
        count = count+1;
    end
end

plot(x,d4, 'k')
M = max(max(M,d4));

xlabel('log(\rho)')
ylabel('Edit Distance')
%title('N=10')
%legend('Low to Medium', 'Low to High', 'Medium to High');
% axis([0 1 0 M+1])
% xticks(0:0.1:1);
% xticklabels({'0','','0.2','','0.4','','0.6','','0.8','','1'})
hold off



