function A=adjMat(particles, BL)
    %returns the adjacency matrix of cluster stored in particles. Bonds
    %exist if distance is <=BL. set to 1 for normal, other for tests. 
    
    %create the particles distance matrix
    D = pDists(particles);
    
    %Initialize adjacency matrix and compute
    tol = 1e-5; 
    L = length(D); 
    A = zeros(L); 
    for i=1:L
        for j=i+1:L
            if D(i,j) < BL + tol
                A(i,j)=1; A(j,i)=1;
            end
        end
    end
end