function  [dist] = comparePartitions(n, a, b, type)
%Creates partitions of merged clusters in the rho descent graphs.

Eout=["LOW","MED","HIGH"];  %List of low, med, high energies
test=[a,b];         %1, 2, or 3 correspond to low, med, high energy test

%Set up total num of clusters for normalization
num_clusters(6)=2; num_clusters(7)=5; num_clusters(8)=13; num_clusters(9)=52;
num_clusters(10)=263; num_clusters(11)=1659;

%Load in the partitions and TF vectors - Morse to Morse
if type == "MM"
    load(strcat('n',num2str(n),Eout(test(1)),'mergeTF.mat'))
    v1=eval(strcat('n',num2str(n),Eout(test(1)),'mergeTF'));
    load(strcat('n',num2str(n),Eout(test(2)),'mergeTF.mat'))
    v2=eval(strcat('n',num2str(n),Eout(test(2)),'mergeTF'));
    load(strcat('n',num2str(n),Eout(test(1)),'partitions.mat'))
    p1=eval(strcat('n',num2str(n),Eout(test(1)),'partitions'));
    load(strcat('n',num2str(n),Eout(test(2)),'partitions.mat'))
    p2=eval(strcat('n',num2str(n),Eout(test(2)),'partitions'));
end

%Load in the partitions and TF vectors - LJ to LJ
if type == "LL"
    load(strcat('n',num2str(n),Eout(test(1)),'LJmergeTF.mat'))
    v1=eval(strcat('n',num2str(n),Eout(test(1)),'LJmergeTF'));
    load(strcat('n',num2str(n),Eout(test(2)),'LJmergeTF.mat'))
    v2=eval(strcat('n',num2str(n),Eout(test(2)),'LJmergeTF'));
    load(strcat('n',num2str(n),Eout(test(1)),'LJpartitions.mat'))
    p1=eval(strcat('n',num2str(n),Eout(test(1)),'LJpartitions'));
    load(strcat('n',num2str(n),Eout(test(2)),'LJpartitions.mat'))
    p2=eval(strcat('n',num2str(n),Eout(test(2)),'LJpartitions'));
end

%Load in the partitions and TF vectors - Morse to LJ
if type == "ML"
    load(strcat('n',num2str(n),Eout(test(1)),'mergeTF.mat'))
    v1=eval(strcat('n',num2str(n),Eout(test(1)),'mergeTF'));
    load(strcat('n',num2str(n),Eout(test(2)),'LJmergeTF.mat'))
    v2=eval(strcat('n',num2str(n),Eout(test(2)),'LJmergeTF'));
    load(strcat('n',num2str(n),Eout(test(1)),'partitions.mat'))
    p1=eval(strcat('n',num2str(n),Eout(test(1)),'partitions'));
    load(strcat('n',num2str(n),Eout(test(2)),'LJpartitions.mat'))
    p2=eval(strcat('n',num2str(n),Eout(test(2)),'LJpartitions'));
end

%Make random trees with same structure as desired tree
if type == "R"
    load(strcat('n',num2str(n),Eout(test(1)),'mergeTF.mat'))
    v1=eval(strcat('n',num2str(n),Eout(test(1)),'mergeTF'));
    load(strcat('n',num2str(n),Eout(test(2)),'LJmergeTF.mat'))
    v2=eval(strcat('n',num2str(n),Eout(test(2)),'LJmergeTF'));
    load(strcat('n',num2str(n),Eout(test(1)),'partitions.mat'))
    p1=eval(strcat('n',num2str(n),Eout(test(1)),'partitions'));
    load(strcat('n',num2str(n),Eout(test(2)),'LJpartitions.mat'))
    p2=eval(strcat('n',num2str(n),Eout(test(2)),'LJpartitions'));
    rho = find(v1==1); rho=rho(end:-1:1);     
    p1= makeRand(p1,n,num_clusters,rho); p2=makeRand(p2,n,num_clusters,rho);
end

%p1{2}
%p2{2}
%Compare merge vectors. If same, perform comparison only at levels with 1
dist=zeros(1,50);  %Vector of distance of partitions at each rho
rho = v1+v2>0; rho = find(rho==1); %Give rho values were either test has merges

%Compute edit distance - symmetrized (divide by 2?)
for i = rho
    a = p1{i};
    b = p2{i};
    dist(i) = compareP(b,a) + compareP(a,b);
end

%fix zero rho entries
for i = 48:-1:2
    if dist(i) == 0
        dist(i) = dist(i+1);
    end
end
end


function s1 = makeRand(p1,n,num_clusters,rho)
%random graph with the same number and length of partitions at each rho
%level

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

function s1 = makeRand2(p1,n,num_clusters,rho)
%random graph with same num of partitions but different lengths, equally
%distibuted
s1={};
for j=rho
    p=p1{j}; %first real partition
    s={}; count=1;
    r = randperm(num_clusters(n));
    bars = sort(unique(randi(num_clusters(n), length(p)-1, 1))); 
    bars = [1; bars; num_clusters(n)+1];
    lengths = circshift(bars, -1) - bars; lengths=lengths(1:end-1);
    for i = 1:length(lengths)
        L=lengths(i);
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
