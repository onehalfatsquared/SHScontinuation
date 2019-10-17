function D=sortedPD(clust)
    %returns a vector of particle distances
    particles=c2p(clust); l=length(particles); D=[]; 
    for i=1:l
        for j=i+1:l
            D=[D,euDist(particles,i,j)];
        end
    end
    D=sort(D);
end