%Creates various functions and distributions related to rhoDescent graph.
%Computes number of clusters as fn of rho, number of merge events as fn of
%rho, and distribution of how many nodes merge in merge event. 
clear

n=10;            %Number of particles   
Eout=["LOW","MED","HIGH"];  %List of low, med, high energies
test=1;         %1, 2, or 3 correspond to low, med, high energy test
Nc=ones(1,50);  %Store number of clusters as fn of rho
Nm=zeros(1,50); %store number of merge events as fn of rho
Nn=cell(1,50);  %store distribution of num of nodes per merge event

%Loop that reads data for each rho and determines merging
for i=49:-1:1
    %Read data from file
    filename=strcat('n',num2str(n),'m',num2str(i),'Sticky',Eout(test),'.txt'); %File name nxrhoy.txt
    fileID=fopen(filename,'r');             %Open file
    formatSpec='%f';                        %Format as float
    data=fscanf(fileID,formatSpec);         %read in data
    fclose(fileID);                         %Close file
    
    %Seperate column vector into an array of clusters, num clusters by 3n.
    num_clusters=length(data)/(3*n);
    clusters=reshape(data,[3*n,num_clusters])';
    
    %First time through, set initial num of clusters and initialize vars
    if i==49
        %Initialize list of merged nodes set to 1
        mergeYet=ones(1,num_clusters);
        
        %Set N_c(50) to be total number of clusters
        Nc(50)=num_clusters; 
    end
    
    %Scale down rho. Test if clusters are same. If so, create a merged
    %node indexed by min cluster, remove duplicates.
    subtractor=0;  %add one every time a node merges
    events=0;      %store num of merge events per rho
    mergeSizes=[]; %store distribution of merge sizes per rho
    for j=1:num_clusters
        kFlag=0; merges=0; 
        if mergeYet(j)==1
            for k=j+1:num_clusters
                %test if j and k are same, if they are non-merged
                if mergeYet(k)==1
                    b=testSame(clusters(j,:),clusters(k,:));
                    if b==1
                        subtractor=subtractor+1; %add one to num of merged nodes
                        merges=merges+1; %add one to num of merges for this j
                        mergeYet(k)=0;  %remove duplicate nodes from list
                        kFlag=1;  %flag to compute num of merge events
                    end
                end
            end
            if kFlag==1
                events=events+1; %add one to num of events
                mergeSizes=[mergeSizes merges]; %add num merges to distr
            end
        end
    end
    
    %store values for rho=i. 
    Nc(i)=Nc(i+1)-subtractor;
    Nm(i)=events;
    Nn{i}=mergeSizes;
    i
end


