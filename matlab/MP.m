function U=MP(particles,rho,E)
    %Evalaute morse potential
    n=length(particles);
    U=0;                     %Initialize sum
    for i=1:n-1
        for j=i+1:n
            Y=exp(-rho*(euDist(particles,i,j)-1));
            U=U+E*Y^2-2*E*Y;
        end
    end
end