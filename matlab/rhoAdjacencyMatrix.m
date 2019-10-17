%Create an adjacency matrix for the rho descent graph. 
clear

n=6;               %Number of particles   
merge=1;           %count number of rho values merged at

%Loop that reads data for each rho and determines merging
for i=49:-1:1
    %Read data from file
    filename=strcat('n',num2str(n),'rho',num2str(i),'TEST.txt'); %File name nxrhoy.txt
    fileID=fopen(filename,'r');             %Open file
    formatSpec='%f';                        %Format as float
    data=fscanf(fileID,formatSpec);      %read in data
    fclose(fileID);                         %Close file
    
    %Seperate column vector into an array of clusters, num clusters by 3n.
    num_clusters=length(data)/(3*n);
    clusters=reshape(data,[3*n,num_clusters])';
    
    if i==49
        num_nodes=num_clusters;        %number of nodes in graph to be updated
        mergeYet=ones(1,num_clusters); %one if not merged, 0 if merged
        nodeIndex=1:num_nodes;         %Gives current index of merged nodes
        A=[];                          %initialize adjacency matrix
    end
    
    %Scale down rho. Test if clusters are same. If so, update the adjadency
    %matrix. 
    for j=1:num_clusters
        mergeList=j;
        if mergeYet(j)==1
            min_clust=j;
            for k=j+1:num_clusters
                %test if j and k are same, if they are non-merged
                if mergeYet(k)==1
                    b=testSame(clusters(j,:),clusters(k,:));
                    if b==1
                        mergeList=[mergeList, k];
                        mergeYet(k)=0;
                    end
                end
            end
            if length(mergeList)>1
                %At least 2 clusters are the same. Merge.
                for nodes=1:length(mergeList)
                    A(nodeIndex(mergeList(nodes)),num_nodes+1)=1;
                    merge=merge+1;
                end
                num_nodes=num_nodes+1;
                nodeIndex(mergeList(1))=num_nodes;
            end
        end
    end
end

A(num_nodes,1)=0; %Make matrix square with zero row at bottom. 


