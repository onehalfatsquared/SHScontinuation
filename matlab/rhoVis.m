%Visualize the merging of clusters as a function of the parameter rho to
%the Morse potential
clear
clf

rng(12,'twister');           %Seed for random color plots
n=6;               %Number of particles   
ticks=50;          %Store the merge points

%Set up for plotting, reverse x axis
hold on
set(gca,'XDir','Reverse')

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
    
    %First time through, compute distances to 1, sort. Create test set.
    %Create axes for plot. x from 50 to 0.5, y from 0.9 to num_clust
    if i==49
        D=zeros(1,num_clusters);
        testSet=1:num_clusters;
        axis([0.5,50,0.9,num_clusters+0.1]);
        for j=2:num_clusters
            [~,~,~,~,~,d]=testSame(clusters(1,:),clusters(j,:),1);
            D(j)=d;
        end
        [D,I]=sort(D);
    end
    
    for j=testSet
        for k=testSet(find(testSet==j)+1:length(testSet))
            %test if j and k are same
            b=testSame(clusters(j,:),clusters(k,:));
            if b==1
                %mark tick at rho, plot traj. up to merging
                ticks=[ticks i];
                y=max(find(I==j),find(I==k));
                clust=I(y);
                y2=min(find(I==j),find(I==k));
                c=rand(1,3);
                line([ticks(1),ticks(length(ticks))],[y,y],'color',c);
                line([ticks(length(ticks)),ticks(length(ticks))],[y,y2],'color',c);
                testSet=testSet(testSet~=clust);
            end
        end
    end
end
hold off

%Plot a horizontal line for the first cluster
line([ticks(1),1],[1,1]);

%Axis titles
xlabel('\rho'); ylabel('Cluster Number');
title(['Cluster Merging Diagram, n=',num2str(n)])

%Plot the tick markers on x giving merge location and on y giving cluster
%number
ticks=unique([1 ticks]); xticks(ticks);
yticks(1:num_clusters); ylbl=I(1:num_clusters); yticklabels(string(ylbl));



