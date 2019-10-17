%get average and std of edit distance of random graphs
clear
a=zeros(1,10);

%num clusters as fn of N
num_clusters(6)=2; num_clusters(7)=5; num_clusters(8)=13; num_clusters(9)=52;
num_clusters(10)=263; num_clusters(11)=1659;

%test params
n = 10;
N = num_clusters(n); %normalizing factor

%do 10 tests
for i=1:10
    a(i) = sum(comparePartitions(n,1,2,"R"))/N;
end

%output
mean(a)
std(a)
