function D=pDists(particles)
    %returns a vector of particle distances
    l=length(particles); 
    for i=1:l
        for j=i+1:l
            D(i,j)=euDist(particles,i,j);
        end
    end
end