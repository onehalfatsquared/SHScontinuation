%plot curves of cluster energy as fn of rho to see hgow they evolve away
%from shs limit

clear;

n = 6; %num particles

%Read my data from file
filename=strcat('n',num2str(n),'adjust.txt');  %File name nxMorse30.txt
fileID=fopen(filename,'r');                     %Open file
formatSpec='%f';                                %Format as float
mydata=fscanf(fileID,formatSpec);               %read in data    
fclose(fileID);                                 %Close file

%Seperate column vector into an array of clusters, num clusters by 3n.
num_clusters=length(mydata)/(3*n);
clusters=reshape(mydata,[3*n,num_clusters])';

%get energies for all clusters in shs limit
bonds = 3*n-6; e0 = -bonds;

%define energy storage
Energies = zeros(num_clusters,50);

%fill column 50 with shs values
Energies(:,50) = e0*ones(num_clusters,1);

%loop over all rho, all clusters, store energy
for rho=1:49
    for i = 1:num_clusters
        c = getCluster(n,rho,i);
        e = MP(c2p(c),rho,1);
        
        Energies(i,rho) = e;
    end
end

% Make figures look better
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',24)

%plot 
figure(1)
hold on
colormap(lines(num_clusters+1));
cmap=colormap;
for i=1:num_clusters
    Plot_color=cmap(i,:);
    plot(log(1:50),Energies(i,:), 'Color', Plot_color)
end
xlabel("log(\rho)");
ylabel("Potential Energy")
%axis([0,50,-30,-10])


















function b = countBonds(cluster)
A = adjMat(c2p(cluster), 1.03);
b = sum(sum(A))/2;
end

