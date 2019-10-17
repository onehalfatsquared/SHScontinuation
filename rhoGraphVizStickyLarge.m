%Create a document to pass to graphviz to visualize merging with rho.
%This version tries to merge at endpoint instead of intermeidate pt
%and uses the sticky files as input. For large files, suppress labels. 
clear

n=8;               %Number of particles   
merge=1;           %count number of rho values merged at
flag=0;            %For seperating rho values in graph
Eout=["LOW","MED","HIGH"];  %List of low, med, high energies
test=2; %1, 2, or 3 correspond to low, med, high energy test

%parameters to the graph
width=0.1;
arrowsize=0.4;
nodesep=0.03;
ranksep=0.3;

%Loop that reads data for each rho and determines merging
for i=49:-1:1
    %Read data from file
    filename=strcat('n',num2str(n),'rho',num2str(i),'Sticky',Eout(test),'.txt'); %File name nxrhoy.txt
    fileID=fopen(filename,'r');             %Open file
    formatSpec='%f';                        %Format as float
    data=fscanf(fileID,formatSpec);      %read in data
    fclose(fileID);                         %Close file
    
    %Seperate column vector into an array of clusters, num clusters by 3n.
    num_clusters=length(data)/(3*n);
    clusters=reshape(data,[3*n,num_clusters])';
    
    %Checkpointing
    i
    
    
    %First time through, make boxes for all clusters with rho=50. 
    %Write output file. Set rank of all rho=50 to be same
    if i==49
        %Store the last merged rho value of each cluster and if merged
        mergeYet=ones(1,num_clusters);
        mergeIndex=50*mergeYet;
        
        %Makr circles for all starting clusters and label. Same rank
        filename=strcat('n',num2str(n),'rhoGraphSticky',Eout(test),'Large.txt'); %File name nxrhoGraph.txt
        fileID=fopen(filename,'w');                     %Open file
        fprintf(fileID, 'digraph Merging_Diagram_n%d  {\n nodesep=%.2f; ranksep=%.2f;\n',n, nodesep, ranksep);
        for j=1:num_clusters
            fprintf(fileID,'"%dp%d" [label="", shape=point  , width=%.2f, regular=1,style=filled] ;\n',j,mergeIndex(j), width);      %Write cluster
        end
        fprintf(fileID,'{rank = same;');
        for j=1:num_clusters
            fprintf(fileID,'"%dp%d";',j,mergeIndex(j));
        end
        fprintf(fileID,'}\n');
        fclose(fileID);
    end
    
    %Scale down rho. Test if clusters are same. If so, create a merged
    %node indexed by min cluster, remove duplicates, set same rank per rho
    rankList=[];
    for j=1:num_clusters
        mergeList=j;
        if mergeYet(j)==1
            min_clust=j; 
            for k=j+1:num_clusters
                %test if j and k are same, if they are non-merged
                if mergeYet(k)==1
                    b=testSame(clusters(j,:),clusters(k,:));
                    if b==1
                        rankList=[rankList min_clust];
                        mergeList=[mergeList, k];
                        mergeYet(k)=0;
                    end
                end
            end
            if length(mergeList)>1
                %At least 2 clusters are the same. Merge.
                filename=strcat('n',num2str(n),'rhoGraphSticky',Eout(test),'Large.txt'); %File name nxrhoGraph.txt
                fileID=fopen(filename,'a+');                     %Open file
                fprintf(fileID,'"%dp%d" [label="", shape=point  , width=%.2f, regular=1,style=filled] ;\n',j,i, width);
                fprintf(fileID, '{"%dp%d"',j,mergeIndex(j));
                for nodes=mergeList(2:end)
                    fprintf(fileID, ', "%dp%d"',nodes,mergeIndex(nodes));
                end
                mergeIndex(j)=i;
                fprintf(fileID, '} -> "%dp%d"[arrowsize=%.2f];\n', j,mergeIndex(j),arrowsize);
                
                fclose(fileID);
                
                
                merge=merge+1;
            end
        end
    end
    if isempty(rankList)==0
        rankList=unique(rankList);
        
        filename=strcat('n',num2str(n),'rhoGraphSticky',Eout(test),'Large.txt'); %File name nxrhoGraph.txt
        fileID=fopen(filename,'a+');                     %Open file
        
        
        
        if flag==1
            fprintf(fileID, '"%dp%d" -> "%dp%d" [style=invis]\n', prev(1), prev(2), rankList(1), mergeIndex(rankList(1)));
        end
        fprintf(fileID,'{rank = same;');
        for j=rankList
            fprintf(fileID,'"%dp%d";',j,mergeIndex(j));
        end
        fprintf(fileID,'}\n');
        fclose(fileID);
        prev=[rankList(1), mergeIndex(rankList(1))];
        flag=1; %Flag to tell if one row has been outputted yet. For organization of layers. 
    end
end

filename=strcat('n',num2str(n),'rhoGraphSticky',Eout(test),'Large.txt'); %File name nxrhoGraph.txt
fileID=fopen(filename,'a+');                     %Open file
fprintf(fileID, '}');                           %end the graphviz file with }
fclose(fileID);                                 %Close file
