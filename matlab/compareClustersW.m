%See if the clusters found by gradient descent / conjugate gradient
%for the Morse potential are included in the list by David Wales. 
clear

n=12;  %Number of particles
rho=30; %rho values
Eout=["LOW","MED","HIGH"];
test=3;

%Read my data from file
filename=strcat('n',num2str(n),'rho',num2str(rho),'Sticky',Eout(test),'.txt');  %File name nxMorse30.txt

fileID=fopen(filename,'r');                     %Open file
formatSpec='%f';                                %Format as float
mydata=fscanf(fileID,formatSpec);               %read in data    
fclose(fileID);                                 %Close file

%Read Wales data from file
filename=strcat('n',num2str(n),'rho',num2str(rho),'Wales.txt');  %File name nxMorse30.txt
fileID=fopen(filename,'r');                    %Open file
formatSpec='%f';                               %Format as float
Wdata=fscanf(fileID,formatSpec);               %read in data    
fclose(fileID);                                %Close file

%Seperate column vector into an array of clusters, num clusters by 3n.
%My clusters
num_clustersM=length(mydata)/(3*n);
clustersM=reshape(mydata,[3*n,num_clustersM])';
%Wales Clusters
num_clustersW=length(Wdata)/(3*n);
clustersW=reshape(Wdata,[3*n,num_clustersW])';

%Test if each of my clusters is contained in Wales clusters
% p=0;   %Plot if 1, dont if 0
% for i=1:num_clustersM
%     count=0; %Check if duplicates
%     clusterM=clustersM(i,:);             %My cluster interested in
%     for j=1:num_clustersW
%        clusterW=clustersW(j,:);         %Wales cluster
%        %Check if pM and pW are the same. if so plot. 
%        %[check,M,W,rmsd]=testSame(clusterM,clusterW);
%        check=testSameHeur2(clusterM,clusterW); 
%        if check==1
%            fprintf('ClusterM %d and ClusterW %d have distance %e\n',i,j,1-check)
%            count=count+1;
%            if p==1
%                M=reshape(M',[3*n,1]); W=reshape(W',[3*n,1]);
%                opts = struct('srad',0.1,'fig',1,'lcolr',0.7*[1,1,1],'lrad',0.03,...
%                    'scolr',1,'salph',0.7,'ifgrid',0,'iftext',1,...
%                    'iftight',1,'az',20,'el',18,'zth',-0.5,'pos',[7,5,15,15]);
%                plotcluster2b(c1,opts);
%                opts = struct('srad',0.1,'fig',2,'lcolr',0.7*[1,1,1],'lrad',0.03,...
%                    'scolr',1,'salph',0.7,'ifgrid',0,'iftext',1,...
%                    'iftight',1,'az',20,'el',18,'zth',-0.5,'pos',[25,5,15,15]);
%                plotcluster2b(c2,opts);
%                pause()
%            end
%            break
%        end
%     end
%     if count>1
%         fprintf('verify\n\n')
%     end
%     if check==0
%         fprintf('ClusterM %d does not match any Wales\n',i)
%     end
% end

%Test if each of Wale's clusters is contained in mine
p=0;   %Plot if 1, dont if 0 
noMatch=[];
for i=1:num_clustersW
    clusterW=clustersW(i,:);             %My cluster interested in
    for j=1:num_clustersM
        if max(abs(clusterW))>2
            check=1;
            fprintf('ClusterW %d is in-admissable\n',i)
            break
        end
        clusterM=clustersM(j,:);         %Wales cluster
       %Check if pM and pW are the same. if so plot. 
       %[check,M,W,rmsd]=testSame(clusterM,clusterW);
       check=testSameHeur2(clusterM,clusterW);
       if check==1
           fprintf('ClusterW %d and ClusterM %d have distance %e\n',i,j,0)
           if p==1
               M=reshape(M',[3*n,1]); W=reshape(W',[3*n,1]);
               opts = struct('srad',0.1,'fig',1,'lcolr',0.7*[1,1,1],'lrad',0.03,...
                   'scolr',1,'salph',0.7,'ifgrid',0,'iftext',1,...
                   'iftight',1,'az',20,'el',18,'zth',-0.5,'pos',[7,5,15,15]);
               plotcluster2b(M,opts);
               opts = struct('srad',0.1,'fig',2,'lcolr',0.7*[1,1,1],'lrad',0.03,...
                   'scolr',1,'salph',0.7,'ifgrid',0,'iftext',1,...
                   'iftight',1,'az',20,'el',18,'zth',-0.5,'pos',[25,5,15,15]);
               plotcluster2b(W,opts);
               pause()
           end
           break
       end
    end
    if check==0
        fprintf('ClusterW %d does not match any of mine\n',i)
        vec=[i; MP(c2p(clusterW),rho,1)];
        noMatch=[noMatch vec];
    end
end

%Test if each of my clusters is unique
% p=0;   %Plot if 1, dont if 0
% unique=0;
% for i=1:num_clustersM
%     cluster1=clustersM(i,:);             %My cluster interested in
%     for j=i+1:num_clustersM
%         cluster2=clustersM(j,:);         %Wales cluster
%        %Check if pM and pW are the same. if so plot. 
%        %[check,M,W,rmsd]=testSame(cluster1,cluster2);
%        check = testSameHeur2(cluster1,cluster2); rmsd=0;
%        if check==1
%            fprintf('ClusterW %d and ClusterW %d have distance %e\n',i,j,rmsd)
%            if p==1
%                M=reshape(M',[3*n,1]); W=reshape(W',[3*n,1]);
%                opts = struct('srad',0.1,'fig',1,'lcolr',0.7*[1,1,1],'lrad',0.03,...
%                    'scolr',1,'salph',0.7,'ifgrid',0,'iftext',1,...
%                    'iftight',1,'az',20,'el',18,'zth',-0.5,'pos',[7,5,15,15]);
%                plotcluster2b(M,opts);
%                opts = struct('srad',0.1,'fig',2,'lcolr',0.7*[1,1,1],'lrad',0.03,...
%                    'scolr',1,'salph',0.7,'ifgrid',0,'iftext',1,...
%                    'iftight',1,'az',20,'el',18,'zth',-0.5,'pos',[25,5,15,15]);
%                plotcluster2b(W,opts);
%                pause()
%            end
%            break
%        end
%     end
%     if check==0
%         unique=unique+1;
%         fprintf('ClusterW %d is unique\n',i)
%     end
% end

%Test if each of Wale's clusters is unique
% p=0;   %Plot if 1, dont if 0
% unique=0;
% for i=1:num_clustersW
%     cluster1=clustersW(i,:);             %My cluster interested in
%     for j=i+1:num_clustersW
%         cluster2=clustersW(j,:);         %Wales cluster
%        %Check if pM and pW are the same. if so plot. 
%        %[check,M,W,rmsd]=testSame(cluster1,cluster2);
%        [check]=testSameHeur2(cluster1,cluster2);
%        if check==1
%            fprintf('ClusterW %d and ClusterW %d have distance %e\n',i,j,1-check)
%            if p==1
%                M=reshape(M',[3*n,1]); W=reshape(W',[3*n,1]);
%                opts = struct('srad',0.1,'fig',1,'lcolr',0.7*[1,1,1],'lrad',0.03,...
%                    'scolr',1,'salph',0.7,'ifgrid',0,'iftext',1,...
%                    'iftight',1,'az',20,'el',18,'zth',-0.5,'pos',[7,5,15,15]);
%                plotcluster2b(M,opts);
%                opts = struct('srad',0.1,'fig',2,'lcolr',0.7*[1,1,1],'lrad',0.03,...
%                    'scolr',1,'salph',0.7,'ifgrid',0,'iftext',1,...
%                    'iftight',1,'az',20,'el',18,'zth',-0.5,'pos',[25,5,15,15]);
%                plotcluster2b(W,opts);
%                pause()
%            end
%            break
%        end
%     end
%     if max(abs(cluster1))>1.7
%         fprintf('ClusterW %d is in-admissable\n',i)
%         check=1;
%     end
%     if check==0
%         unique=unique+1;
%         fprintf('ClusterW %d is unique\n',i)
%     end
% end
