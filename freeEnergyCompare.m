%Evaluate the free energy of all the unique clusters for a given number of
%particles, range, sticky parameter, and temperature
clear

%parameters
n=7;  %Number of particles
rho=49; %range parameter
kap=20; %medium sticky parameter value
E=stickyNewton(8.3, rho, kap); %find E to make kappa constant
T=5; %temp in kelvin

%Read my data from file
filename=strcat('n',num2str(n),'rho',num2str(rho),'StickyMED.txt');  %File name nxMorse30.txt
fileID=fopen(filename,'r');                     %Open file
formatSpec='%f';                                %Format as float
mydata=fscanf(fileID,formatSpec);               %read in data    
fclose(fileID);                                 %Close file

%Seperate column vector into an array of clusters, num clusters by 3n.
num_clusters=length(mydata)/(3*n);
clusters=reshape(mydata,[3*n,num_clusters])';

%find the unique clusters for given value of rho
num_unique=0; unique=[]; 
for i=1:num_clusters
    cluster1=clusters(i,:);             
    for j=i+1:num_clusters
        cluster2=clusters(j,:);         
       %Check if clusters are the same  
       [check,~,~,rmsd]=testSame(cluster1,cluster2);
       if check==1
           fprintf('Cluster %d and Cluster %d have distance %e\n',i,j,rmsd)
           break
       end
    end
    if check==0
        %add unique cluster
        num_unique=num_unique+1;
        fprintf('Cluster %d is unique\n',i)
        unique=[unique i];
    end
end

%if all clusters are the same then add one to unique list. 
if isempty(unique)
    unique=1; 
end

F=freeEnergyM([1,2,3,4,5],clusters(unique,:),rho,E,T)
%F=freeEnergyM([3,5],clusters(unique,:),rho,E,T)
%F=freeEnergyM([4],clusters(unique,:),rho,E,T)




