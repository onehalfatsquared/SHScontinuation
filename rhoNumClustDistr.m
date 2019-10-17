%Computes the number of unique clusters for fixed n, as a function of rho. 
clear

n=7;               %Number of particles   
Eout=["LOW","MED","HIGH"];  %List of low, med, high energies
test=1; %1, 2, or 3 correspond to low, med, high energy test
Nc=ones(1,50);    %Initialize storage for N_c(rho)

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
    
    %First time through, set initial num of clusters and initialize vars
    if i==49
        %Initialize list of merged nodes set to 1
        mergeYet=ones(1,num_clusters);
        
        %Set N_c(50) to be total number of clusters
        Nc(50)=num_clusters;
    end
    
    %Scale down rho. Test if clusters are same. If so, create a merged
    %node indexed by min cluster, remove duplicates.
    subtractor=0; %Increment to get num nodes to subtract when merge
    for j=1:num_clusters
        if mergeYet(j)==1
            for k=j+1:num_clusters
                %test if j and k are same, if they are non-merged
                if mergeYet(k)==1
                    b=testSame(clusters(j,:),clusters(k,:));
                    if b==1
                        mergeYet(k)=0; %remove duplicates
                        subtractor=subtractor+1;
                    end
                end
            end
        end
    end
    Nc(i)=Nc(i+1)-subtractor;
end

