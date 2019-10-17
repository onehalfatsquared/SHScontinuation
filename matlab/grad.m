function g=grad(cluster,rho,E)
    %Compute the gradient of a morse potential at initial point cluster
    g=cluster;                      %Initialize gradient vector 
    c=length(g);                    %Number of coordinates
    particles=c2p(cluster);         %Store as particle array
    %Compute the gradient. Ignore coordinates 1,2,3,5,6,9. 
    for i=4:c
        if i~=5 && i~=6 && i~=9
            dU=0;                   %Initialize sum
            j=ceil(i/3);            %particle containing coordinate i
            elements=1:c/3; elements(j)=[]; %Sum over all particles neq j
            for k=elements
                Y=exp(-rho*(euDist(particles,j,k)-1));
                dU=dU+(Y-Y^2)*(cluster(i)-cluster(i-3*(j-k)))/euDist(particles,j,k);
            end
            g(i)=dU*2*E*rho;
        end
    end
end