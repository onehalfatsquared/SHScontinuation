function U=LJ(particles,m,E)
    %Evalaute Lenard Jones potential with n=2m 
    n=length(particles);
    U=0;                     %Initialize sum
    for i=1:n-1
        for j=i+1:n
            R=euDist(particles,i,j)^(-m);
            U=U+E*(R^2-2*R);
        end
    end
end