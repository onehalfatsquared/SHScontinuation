%find the first time a cluster merges and store it
clear

n=11;               %Number of particles   
merge=1;           %count number of rho values merged at
flag=0;            %For seperating rho values in graph
Eout=["LOW","MED","HIGH"];  %List of low, med, high energies
test=2; %1, 2, or 3 correspond to low, med, high energy test

%Loop that reads data for each rho and determines merging
for i=49:-1:1
    i
    %Read data from file
    filename=strcat('n',num2str(n),'rho',num2str(i),'Sticky',Eout(test),'.txt'); %File name nxrhoy.txt
    fileID=fopen(filename,'r');             %Open file
    formatSpec='%f';                        %Format as float
    data=fscanf(fileID,formatSpec);      %read in data
    fclose(fileID);                         %Close file
    
    %Seperate column vector into an array of clusters, num clusters by 3n.
    num_clusters=length(data)/(3*n);
    clusters=reshape(data,[3*n,num_clusters])';
    
    
    %First time through, make boxes for all clusters with rho=50. 
    %Write output file. Set rank of all rho=50 to be same
    if i==49
        %Store the last merged rho value of each cluster and if merged
        mergeYet=ones(1,num_clusters);
        mergeIndex=50*mergeYet;
        
        %create the first merge cell
        for s=1:num_clusters
            fM{s,1}=s; fM{s,2}=0; %cluster and first merge rho. set to 0 if not merged yet.
        end
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
                    b=testSameHeur2(clusters(j,:),clusters(k,:));
                    if b==1
                        rankList=[rankList min_clust];
                        mergeList=[mergeList, k];
                        mergeYet(k)=0;
                        if fM{j,2} == 0
                            fM{j,2} = i;
                        end
                        if fM{k,2} == 0
                            fM{k,2} = i;
                        end
                    end
                end
            end
            if length(mergeList)>1
                %At least 2 clusters are the same. Merge.
                mergeIndex(j)=i;
                merge=merge+1;
            end
        end
    end
    if isempty(rankList)==0
        rankList=unique(rankList);
    end
end

outfile = strcat('n',num2str(n),'fM.mat'); %File name nxrhoy.txt
save(outfile, 'fM'); 




