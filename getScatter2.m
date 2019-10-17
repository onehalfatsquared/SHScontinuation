%get data for jump scatter plot - with rmin
clear

%set parameters
n = 11; %done up to n=9, do not redo

%load the file with rho marge values and jump
%filename = strcat('n',num2str(n),'fM.mat'); %File name nxrhoy.txt
%load(filename)
%load('jump2')

%get number of clusters
%L = length(fM)

%initialize Jump vector
% jump2 = {};
% for i=1:50
%     jump2{i,1} = i;
%     jump2{i,2} = [];
%     jump2{i,3} = []; 
% end

%get the data
% for i =1:L
%     i
%     %get the position of the cluster at rho value just before merge and
%     %after
%     rho = fM{i,2}+1;
%     cb = getCluster(n,rho,i); 
%     ca = getCluster(n,rho-1,i);
%     cr = getCluster(n,49,i);
%     
%     %get rmsd and rmin
%     rm = rmin(cr, rho);
%     [~,~,~,~,~,d] = testSame(cb,ca,1);
%     
%     jump2{rho-1,2} = [jump2{rho-1,2} d];
%     jump2{rho-1,3} = [jump2{rho-1,3} rm]; 
% end
% 
% save('jump2.mat','jump2')
i = 4
c = getCluster(n, 49, i);
d = getCluster(n, 30, i);
%rmin(c, 49)
figure(1)
plotcluster2b(c')
figure(2)
plotcluster2b(d')




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