%Create a  text document to pass to python library that computes edit 
%distances between graphs. 
clear

n=10;               %Number of particles   
flag=0;            %For seperating rho values in graph
Eout=["LOW","MED","HIGH"];  %List of low, med, high energies
test=1; %1, 2, or 3 correspond to low, med, high energy test

%Loop that reads data for each rho and determines merging
for i=49:-1:1
    %Read data from file
    filename=strcat('n',num2str(n),'m',num2str(i),'Sticky',Eout(test),'.txt'); %File name nxrhoy.txt
    fileID=fopen(filename,'r');             %Open file
    formatSpec='%f';                        %Format as float
    data=fscanf(fileID,formatSpec);      %read in data
    fclose(fileID);                         %Close file
    
    %Seperate column vector into an array of clusters, num clusters by 3n.
    num_clusters=length(data)/(3*n);
    clusters=reshape(data,[3*n,num_clusters])';
    
    %Checkpointing
    i
    
    %First time through, create all the child nodes at the top of the graph
    if i==49
        %Store the last merged rho value of each cluster and if merged
        mergeYet=ones(1,num_clusters);
        mergeIndex=50*mergeYet;
        
        %Make child nodes
        filename=strcat('n',num2str(n),'mGraphSticky',Eout(test),'PY.py'); %File name nxrhoGraph.txt
        fileID=fopen(filename,'w');                     %Open file
        fprintf(fileID, 'from zss import Node\n\ndef createTree():\n');
        for j=1:num_clusters
            fprintf(fileID,'\t%s%dp%d = Node("%dp%d", [])\n',Eout(test),j,mergeIndex(j),j,mergeIndex(j));      %Write cluster
        end
        fclose(fileID);
    end
    
    %Scale down rho. Test if clusters are same. If so, create a merged
    %node indexed by min cluster, remove duplicates, set same rank per rho
    for j=1:num_clusters
        mergeList=j;
        if mergeYet(j)==1 
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
                filename=strcat('n',num2str(n),'mGraphSticky',Eout(test),'PY.py'); %File name nxrhoGraph.txt
                fileID=fopen(filename,'a+');                     %Open file
                saveLHS=strcat(Eout(test),num2str(j),'p',num2str(i));
                fprintf(fileID, '\t%s%dp%d = Node("%dp%d", [%s%dp%d', Eout(test),j,i,j,i,Eout(test),j,mergeIndex(j));
                for nodes=mergeList(2:end)
                    fprintf(fileID, ', %s%dp%d',Eout(test),nodes,mergeIndex(nodes));
                end
                fprintf(fileID, '])\n');
                mergeIndex(j)=i;
                fclose(fileID);
            end
        end
    end
end

filename=strcat('n',num2str(n),'mGraphSticky',Eout(test),'PY.py'); %File name nxrhoGraph.txt
fileID=fopen(filename,'a+');               %Open file
fprintf(fileID, '\treturn %s',saveLHS);    %create root
fclose(fileID);                            %Close file