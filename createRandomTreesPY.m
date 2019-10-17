%Create random trees with same number of levels and nodes per level as the
%actual rho descent trees to get an idea of how different the energy
%landscapes are
clear

n=11;               %Number of particles
num_clusters(6)=2; num_clusters(7)=5; num_clusters(8)=13; num_clusters(9)=52; 
num_clusters(10)=263; num_clusters(11)=1659; 

%Create the vertical distribution of nodes
load(strcat('Nmn',num2str(n),'StickyMED.mat'))
distrib=eval(strcat('Nmn',num2str(n),'StickyMED'));
distrib=distrib(distrib~=0); 
distrib=[distrib num_clusters(n) ];
%distrib=[1,1,5]; %for n=7

%Make random graph 1
totalN=sum(distrib); layers=length(distrib); 
filename=strcat('n',num2str(n),'randomGraphsPY.py'); %File name nxrhoGraph.txt
fileID=fopen(filename,'w');                     %Open file
fprintf(fileID, 'from zss import Node\n\ndef createTree1():\n');
for j=1:totalN
    fprintf(fileID,'\tn%d = Node("%d", [])\n',j,j);      %Write cluster
end

for i=2:layers
    num_below= sum(distrib(1:i-1)); %Num of nodes below this layer
    nodes=(1:distrib(i))+num_below; %Index of nodes at this layer
    for j=1:length(nodes)
        r=randi(num_below);
        fprintf(fileID,'\tn%d.addkid(n%d)\n',r,nodes(j));      %Write cluster
    end
end

fprintf(fileID,'\treturn n1') ;     %Write cluster


%Make random graph 2
fprintf(fileID, '\n\n\ndef createTree2():\n');
for j=1:totalN
    fprintf(fileID,'\tn%d = Node("%d", [])\n',j,j);      %Write cluster
end

for i=2:layers
    num_below= sum(distrib(1:i-1)); %Num of nodes below this layer
    nodes=(1:distrib(i))+num_below; %Index of nodes at this layer
    for j=1:length(nodes)
        r=randi(num_below);
        fprintf(fileID,'\tn%d.addkid(n%d)\n',r,nodes(j));      %Write cluster
    end
end

fprintf(fileID,'\treturn n1') ;     %Write cluster
fclose(fileID);