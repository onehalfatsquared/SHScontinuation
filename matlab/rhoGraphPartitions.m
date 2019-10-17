%Creates partitions of merged clusters in the rho descent graphs. 
clear

n=8;            %Number of particles   
Eout=["LOW","MED","HIGH"];  %List of low, med, high energies
test=3;         %1, 2, or 3 correspond to low, med, high energy test
Partitions=cell(1,50);  %store partitions for each rho
mergeTF=zeros(1,50);    %1 if merge at index i, 0 else

%Loop that reads data for each rho and determines merging
for i=49:-1:1
    %Read data from file
    filename=strcat('n',num2str(n),'rho',num2str(i),'Sticky',Eout(test),'.txt'); %File name nxrhoy.txt
    fileID=fopen(filename,'r');             %Open file
    formatSpec='%f';                        %Format as float
    data=fscanf(fileID,formatSpec);         %read in data
    fclose(fileID);                         %Close file
    
    %Seperate column vector into an array of clusters, num clusters by 3n.
    num_clusters=length(data)/(3*n);
    clusters=reshape(data,[3*n,num_clusters])'; 
    
    %First time through, set singleton partition
    if i==49
        for j=1:num_clusters
            p{j}=j;
        end
        Partitions{50}=p; 
        
        %Initialize list of merged nodes set to 1
        mergeYet=ones(1,num_clusters);
    end
    
    %Scale down rho. Test if clusters are same. If so, create a merged
    %node indexed by min cluster, remove duplicates.
    for j=1:num_clusters
        if mergeYet(j)==1
            for k=j+1:num_clusters
                %test if j and k are same, if they are non-merged
                if mergeYet(k)==1
                    b=testSameHeur2(clusters(j,:),clusters(k,:));
                    if b==1
                        p = updatePartition(p,j,k);
                        mergeYet(k)=0;  %remove duplicate nodes from list
                        mergeTF(i)=1;   %there is merge event at this rho
                    end
                end
            end
        end
    end

    %store values for rho=i. 
    Partitions{i}=p;
    i
end

function p = updatePartition(p,j,k)
    %Update partition at merge event
    for i=1:length(p)
        if any(p{i}==k)
            index1=i;
        elseif any(p{i}==j)
            index2=i;
        end
    end
    q=p{index1}; p(index1)=[];
    p{index2}=sort([p{index2}, q]);   
end

