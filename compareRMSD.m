%get rmsd between clusters wqith different kappa / different potential
clear

n=9;
num_clusters(6)=2; num_clusters(7)=5; num_clusters(8)=13; num_clusters(9)=52;
num_clusters(10)=263; num_clusters(11)=1659;

a = getCluster2(n,40,4,"MED", "rho");
b = getCluster2(n,40,4,"LOW", "rho");

vals = cell(49,1);

N = num_clusters(n);
for rho = 1:49
    for c = 1:N
        a = getCluster2(n,rho,c,"HIGH", "rho");
        b = getCluster2(n,rho,c,"HIGH", "m");
        %[~,~,~,~,~,d]=testSame(a,b,1);
        d = rmsd(a,b);
        vals{rho} = [vals{rho},d];
    end
end
    
rhop = [];
mark = [];

for i = 1:length(vals)
    q=vals{i};
    for j=1:length(q)
        rhop = [rhop i];
        mark = [mark q(j)];
    end
    rho_avg(i) = i;
    mark_avg(i) = mean(vals{i});
    e(i) = std(vals{i});
end

figure(1)
scatter(rhop,mark)
xlabel('log(\rho)')
ylabel('rmsd')

figure(2)
scatter(rhop,mark)
xlabel('\rho')
ylabel('rmsd')

axes('Position',[.55 .55 .35 .35])
box on
scatter(log(rhop),mark)
xlabel('log(\rho)')
ylabel('rmsd')

% errorbar(1./rho_avg, mark_avg, e)
% xlabel('1/\rho')
% ylabel('rmsd')
%title('N=9, Morse MED to HIGH')
%title('N=9, MED Morse to LJ');













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