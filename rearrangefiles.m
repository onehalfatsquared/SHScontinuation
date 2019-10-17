%Rearrange data sets such that the first three particles are fixed such
%that coordinates 1,2,3,5,6,9 are always 0. 

clear

%Read data from file
n=12;                                    %Number of particles
filename=strcat('n',num2str(n),'.txt'); %files are nx.txt s.t. x= num part.
fileID=fopen(filename,'r');             %Open file
formatSpec='%f';                        %Format as float
data=fscanf(fileID,formatSpec);         %read in data    
fclose(fileID);                         %Close file

%Seperate column vector into an array of clusters, num clusters by 3n.
num_clusters=length(data)/(3*n);
clusters=reshape(data,[3*n,num_clusters])';

%Rearrange and write to a new file
filename=strcat('n',num2str(n),'adjust.txt'); %files are nx.txt s.t. x= num part.
fileID=fopen(filename,'w');             %Open file
for q=1:num_clusters
    %Seperate a cluster into particles
    c_num=q;                               %Index of cluster interested in
    cluster=clusters(c_num,:);             %Cluster interested in
    particles=zeros(n,3);                  %Store particle locations of cluster
    for i=1:n
        particles(i,:)=clusters(c_num,3*i-2:3*i);
    end
    
    %arrange such that first particle is (0,0,0), second is (x,0,0), and
    %third is (x,x,0);
    a=[0,0,0];
    rowa=find(ismember(particles(:,1:3),a,'rows'),1);
    particles([1,rowa],:)=particles([rowa,1],:);
    b=[0,0];
    rowb=find(ismember(particles(2:n,2:3),b,'rows'),1)+1;
    particles([2,rowb],:)=particles([rowb,2],:);
    c=[0];
    rowc=find(ismember(particles(3:n,3:3),c,'rows'),1)+2;
    particles([3,rowc],:)=particles([rowc,3],:);
    
    %convert array to vector
    new_cluster=[];
    for i=1:n
        new_cluster=[new_cluster particles(i,:)];
    end
    
    %Write to a new file
    fprintf(fileID,'%.16f \n',new_cluster);           %Write cluster
end
fclose(fileID);                         %Close file 