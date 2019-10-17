function [b]=testSameHeur(c1,c2)
    %Test if clusters 1 and 2 are the same up to translation, rotation, and
    %permutation. Uses interparticle dist heuristic  
    
    %Compute interparticle distances of each cluster. Compare to tol
    d1=sortedPD(c1);
    d2=sortedPD(c2);
    distTol=1e-5;
    
    if max(abs(d1-d2))>distTol
        b=0;
    else
        b=1;
    end
end
