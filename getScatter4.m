%get data for jump scatter plot - with rmin and cluster index
clear

%set parameters
n = 10; %done up to n=9, do not redo

%load the file with rho marge values and jump
filename = strcat('n',num2str(n),'fM2.mat'); %File name nxrhoy.txt
load(filename)
load('jumpSP')

%get number of clusters
[L,~] = size(mg)

% initialize Jump vector
%jumpSP{1,1} = [];
% for i=1:50
%     jumpSP{i,1} = [];
% end


%get the data
for i =1:L
    i
    %get the position of the cluster at rho value just before merge and
    %after
    rho = mg{i,1}+1;
    jumps=[];
    for clusts = mg{i,2}
        cb = getCluster(n,rho,clusts); 
        ca = getCluster(n,rho-1,clusts);
        [~,~,~,~,~,d] = testSame(ca,cb,1);
        jumps = [jumps d];
    end
    jumpSP{end+1,1} = [jumps];
    jumpSP{end,2} = rho -1;
    
   
end

save('jumpSP.mat','jumpSP')




function rc = rmin(cluster, rho)
    %get rmin value
    comp_val = [1.01, 1.06];
    if rho>30
        comp = 1;
    else
        comp = 2;
    end
    
    di=unique(pDists(c2p(cluster))); di=di(di>comp_val(comp)); %use 1.01 for rho 38, 1.05 for rho 6
    rc = min(di);
end

function x=rmsd(c1,c2)
    %compute rmsd b/w clusters, assuming they are already in the correct
    %permuatiton
    
    %Convert clusters to particle arrays. 
    p1=c2p(c1); p2=c2p(c2);
    
    %Subtract the "Center of Mass" from each set of particles
    p1=p1-ones(size(p1))*diag(sum(p1)/length(p1));
    p2=p2-ones(size(p2))*diag(sum(p2)/length(p2));
    
    %align the last particles
    u = p1(end,:); v = p2(end,:);
    
    %compute rotation matrix
    R = vecRot(u,v); 
    
    %apply to every particle
    for i=1:length(p1)
        p1(i,:)=p1(i,:)*R'; %Apply rotation to each particle
    end
    
    %optimal rotation
    try
        U=findRot(p1,p2);     %Optimal rotation of p1 to p2
    catch
        x=0;
        return
    end
    
    %compute rmsd
    x = sqrt(1/length(p1)*sum(sum((p1-p2).^2)));    %L2 distance
end