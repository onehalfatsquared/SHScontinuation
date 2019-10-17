function [b]=testSameHeur2(c1,c2)
    %Test if clusters 1 and 2 are the same up to translation, rotation, and
    %permutation. Optionally, determine distance between clusters if arg3 given.  
    
    %parameters
    b=0; %default to not same
    potTol=1e-4; distTol=1e-4; 
    
    %Convert clusters to particle arrays. 
    p1=c2p(c1); p2=c2p(c2);
    
    %First compute energy to see if clusters are potentially the same
    rho=1; E=1;
    E1=MP(p1,rho,E); E2=MP(p2,rho,E);
    if abs(E1-E2)>potTol
        return
    end
    
    %Subtract the "Center of Mass" from each set of particles
    p1=p1-ones(size(p1))*diag(sum(p1)/length(p1));
    p2=p2-ones(size(p2))*diag(sum(p2)/length(p2));
    
    %Compute distances to origin of each particle
    d1=sqrt(sum(p1.^2,2));
    d2=sqrt(sum(p2.^2,2));
    
    %Sort these lists of distances. Compare within tol. If same, permute
    %particle lists
    %[sd1,I1]=sort(d1), [sd2,I2]=sort(d2), max(abs(sd1-sd2)) %use this line to print distance lists and diff
    [sd1]=sort(d1); [sd2]=sort(d2); 
    if max(abs(sd1-sd2))>distTol
        return
    end
    b=1; return   %This line returns b=1 if list of distances is the same
    
end