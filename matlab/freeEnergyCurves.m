%plot curves of cluster free energy as fn of rho to see hgow they evolve away
%from shs limit

clear;

n = 8; %num particles
beta = 1;

%Read my data from file
filename=strcat('n',num2str(n),'adjust.txt');  %File name nxMorse30.txt
fileID=fopen(filename,'r');                     %Open file
formatSpec='%f';                                %Format as float
mydata=fscanf(fileID,formatSpec);               %read in data    
fclose(fileID);                                 %Close file

%Seperate column vector into an array of clusters, num clusters by 3n.
num_clusters=length(mydata)/(3*n);
clusters=reshape(mydata,[3*n,num_clusters])';

%get potential energies for all clusters in shs limit
bonds = 3*n-6; e0 = -bonds;

%define energy storage
Energies = zeros(num_clusters,50);

%define probability storage
p = zeros(1,num_clusters);

%fill column 50 with shs values
for i=1:num_clusters
    clust = clusters(i,:);
    Zr = rotPF(clust);
    Zv = vibPF(clust,50);
    boltz = exp(-beta*e0);
    p(i) = boltz*Zr*Zv;
end
Z = sum(p);
P = p/Z;
Energies(:,50) = -1/beta*log(P);
    
    

%loop over all rho, all clusters, store energy
for rho=1:49
    %get probability of each cluster un-normalized
    for i = 1:num_clusters
        clust = getCluster(n,rho,i);
        Zr = rotPF(clust);
        Zv = vibPF(clust,rho);
        boltz = exp(-beta*e0);
        p(i) = boltz*Zr*Zv;
    end
    
    %look for repeat clusters to throw away in normalization
    num_unique=0; unique=[]; 
    for i=1:num_clusters
        cluster1=getCluster(n,rho,i);            
        for j=i+1:num_clusters
            cluster2=getCluster(n,rho,j);          
           %Check if clusters are the same  
           [check,~,~,rmsd]=testSame(cluster1,cluster2);
           if check==1
               p(j) = 0;
               break
           end
        end
    end

    %normalize, compute FE
    Z = sum(p);
    P = p/Z;
    Energies(:,rho) = -1/beta*log(P);
    
    
end

% Make figures look better
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',24)

%plot 
figure(1)
hold on
for i=1:num_clusters
    plot(1:50,Energies(i,:))
end
xlabel("\rho");
ylabel("Free Energy")
axis([0,50,-1,10])




function P = computeProb(clust,rho)
    %Compute probability up to constants of a cluster
    
    beta = 1;
    Zr = rotPF(clust);
    Zv = vibPF(clust,rho,1);
    Z=Zr*Zv; %partition function up to constants
    U=MP(c2p(clust),100,1); %morse potential energy
    boltz=exp(-beta*U); %boltzmann factor
    P=Z*boltz; %probability, un-normalized
end

function zR = rotPF(clust)
    %evaluate the rotational partition function of a cluster
    
    M=inertiaTensor(clust); I=det(M); %inertia tensor + determinant
    s=symNum(clust); %evaluate symmetry number
    zR=sqrt(I)/s;  %eval partition fn
end

function zV = vibPF(clust,rho)
    %evaluate vibrational partition function
    
    n=length(clust)-6; %num of Dofs
    H=hessMorse(clust,rho,1);  %hessian with zeros eigs
    e=eig(H); e=e(e>0); e=prod(e); %remove zero eigs, take product
    zV=(2*pi)^(n/2)/sqrt(e); %evaluate vib partition fn
end

function s = symNum(clust)
    %evaluate the symmetry number of a cluster
    
    tol=1e-5; %tolerance for same matrix elements
    particles=c2p(clust); parts=length(particles); % particle array and num particles
    ipdMatrix = pDists(particles); % inter particle distance matrix
    AP=perms(1:parts); %list of all perms
    s=0; %initialize sym number
    %Loop to check how many perms give same ipd matrix
    for i=1:factorial(parts)
        test=particles(AP(i,:),:);
        testMatrix=pDists(test); 
        if max(max(abs(testMatrix-ipdMatrix)))<tol
            s=s+1;
        end
    end
end


















function b = countBonds(cluster)
A = adjMat(c2p(cluster), 1.03);
b = sum(sum(A))/2;
end

