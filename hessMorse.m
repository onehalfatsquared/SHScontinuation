function H=hessMorse(cluster,rho,E)
    %Compute the Hessian of the Morse rho potential for given cluster
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
                        Y=exp(-rho*(r-1));
                        %Compute term 1
                        T1=(eq(i,s)-eq(i-3*(j-k),s))/r*Y*(1-Y);
                        %Compute term 2
                        if j==ceil(s/3)
                            P2=-1/r^3*(cluster(s)-cluster(s-3*(j-k)));
                        else
                            P2=0;
                            if eq(k,ceil(s/3))
                                P2=-1/r^3*(cluster(s)-cluster(s-3*(k-j)));
                            end
                        end
                        T2=(cluster(i)-cluster(i-3*(j-k)))*Y*(1-Y)*P2;
                        %Compute term 3
                        P3=rho*Y*r^2*P2;
                        T3=(cluster(i)-cluster(i-3*(j-k)))/r*(1-Y)*P3;
                        %Compute term 4
                        T4=(cluster(i)-cluster(i-3*(j-k)))/r*Y*(-P3);
                        d2U=d2U+T1+T2+T3+T4;
                    end
                    H(i,s)=2*E*rho*d2U;
                end
            end
        end
    end
end
