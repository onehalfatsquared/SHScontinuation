function g=LJgrad(cluster,m,E)
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
                r=euDist(particles,j,k);
                dU=dU+(cluster(i)-cluster(i-3*(j-k)))*r^(-m-2)*(1-r^(-m));
            end
            g(i)=dU*2*E*m;
        end
    end
end