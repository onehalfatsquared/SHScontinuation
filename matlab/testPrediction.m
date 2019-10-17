%Script to test if which clusters merging can be predicted.
clear

%Read data from file
n=10;                                    %Number of particles
filename=strcat('n',num2str(n),'adjust.txt'); %files are nx.txt s.t. x= num part.
fileID=fopen(filename,'r');             %Open file
formatSpec='%f';                        %Format as float
data=fscanf(fileID,formatSpec);         %read in data    
fclose(fileID);                         %Close file

%Seperate column vector into an array of clusters, num clusters by 3n.
num_clusters=length(data)/(3*n);
clusters=reshape(data,[3*n,num_clusters])';

%Choose clusters to compare - loop
test = 41;
testc = clusters(test,:);

for i = 50:num_clusters
    if i ~= test
        b = predictMerge(testc, clusters(i,:));
        fprintf('Cluster %d and cluster %d return a %d\n', test, i, b)
    end
end
    

%34 and 58