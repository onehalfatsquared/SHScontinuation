function H=hessLJ(cluster,m,E)
    %Compute the Hessian of the LJ potential for given cluster
    parts=length(cluster);  %Number of particles
    H=zeros(parts);           %Initiliaze Hessian
    particles=c2p(cluster); %Particle array
    %Compute the matrix elements. Ignore coordinates 1,2,3,5,6,9.
    for i=4:parts
        if i~=5 && i~=6 && i~=9
            for s=4:parts
                if s~=5 && s~=6 && s~=9
                    d2U=0;  %Initialize the sum
                    j=ceil(i/3);
                    elements=1:parts/3; elements(j)=[]; %Sum over all particles neq j
                    for k=elements
                        r=euDist(particles,j,k);
                        %Compute term 1
                        T1=(eq(i,s)-eq(i-3*(j-k),s))*r^(-m-2)*(1-r^(-m));
                        %Compute term 2
                        if j==ceil(s/3)
                            P2=-(m+2)*r^(-m-4)*(cluster(s)-cluster(s-3*(j-k)));
                            P3=m*r^(-m-2)*(cluster(s)-cluster(s-3*(j-k)));
                        else
                            P2=0; P3=0;
                            if eq(k,ceil(s/3))
                                P2=-(m+2)*r^(-m-4)*(cluster(s)-cluster(s-3*(k-j)));
                                P3=m*r^(-m-2)*(cluster(s)-cluster(s-3*(k-j)));
                            end
                        end
                        T2=(cluster(i)-cluster(i-3*(j-k)))*P2*(1-r^(-m));
                        %Compute term 3
                        T3=(cluster(i)-cluster(i-3*(j-k)))*r^(-m-2)*P3;
                        d2U=d2U+T1+T2+T3;
                    end
                    H(i,s)=2*E*m*d2U;
                end
            end
        end
    end
end
