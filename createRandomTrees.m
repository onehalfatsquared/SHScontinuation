%Create random trees with same number of levels and nodes per level as the
%actual rho descent trees to get an idea of how different the energy
%landscapes are
clear

n=11;               %Number of particles
Eout = ["LOW","MED","HIGH"];
test=3; 
num_clusters(6)=2; num_clusters(7)=5; num_clusters(8)=13; num_clusters(9)=52; 
num_clusters(10)=263; num_clusters(11)=1659; 

%Create the vertical distribution of nodes
load(strcat('n',num2str(n),Eout(test(1)),'mergeTF.mat'))
v1=eval(strcat('n',num2str(n),Eout(test(1)),'mergeTF'));
load(strcat('n',num2str(n),Eout(test(1)),'partitions.mat'))
p1=eval(strcat('n',num2str(n),Eout(test(1)),'partitions'));

rho = find(v1==1); rho=rho(end:-1:1);

% p=p1{rho(1)}; %first real partition
% r = randperm(num_clusters(n)); 
% s={}; count=1; s1={};
% for i = 1:length(p); 
%     L=length(p{i});
%     s{i}=r(count:count+L-1);
%     count=count+L;
% end
% s1{rho(1)}=s; 
% 
% for i=rho(2:end)
%     p=p1{i};
%     for j=1:length(p)
%         L=length(p({i})); 
%         diff=L-length(s{i});
%         if diff~=0
%             
%         


a= makeRand(p1,n,num_clusters,rho); b=makeRand(p1,n,num_clusters,rho);
dist=zeros(1,50);  %Vector of distance of partitions at each rho
rho = v1>0; rho = find(rho==1); %Give rho values were either test has merges

for i = rho
    Q1 = a{i};
    Q2 = b{i};
    dist(i) = compareP(Q1,Q2);
end

%dist
sum(dist)/num_clusters(n)


function s1 = makeRand(p1,n,num_clusters,rho)
% r = randperm(num_clusters(n));
s1={};
for j=rho
    p=p1{j}; %first real partition
    s={}; count=1;
    r = randperm(num_clusters(n));
    for i = 1:length(p)
        L=length(p{i});
        s{i}=r(count:count+L-1);
        count=count+L;
    end
    s1{j}=s;
end
end

function e = compareP(a,b)
    %Compare the partitions a and b using family distance metric
    e=0;   %initialize edit distance
    for i=1:length(a)
        if length(a{i})>1
            pairs=nchoosek(a{i},2);
            [m,~]=size(pairs);
            for j=1:m
                x=pairs(j,1); y=pairs(j,2); 
                for k=1:length(b)
                    if any(b{k}==x)
                        if any(b{k}==y)
                            break
                        else
                            e=e+1;
                            break
                        end
                    end
                end
            end
        end
    end 
end
