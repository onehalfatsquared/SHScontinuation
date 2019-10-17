function [cluster,counter]=sfNewton(cluster,rho,E)
    %Saddle free newton iterations to find minima

    cutoff=100;   %Cut off if doesn't converge
    counter=0;    %Count number of iterations
    tol=1e-8;    %Tolerance for updates
    
    while counter<cutoff
        H=hessMorse(cluster,rho,E);
        %Delete zero rows and columns 1,2,3,5,6,9.
        H( ~any(H,2), : ) = [];  %rows
        H( :, ~any(H,1) ) = [];  %columns
        
        g0=grad(cluster,rho,E);
        alt_clust=cluster;
        g0(9)=[];g0(6)=[];g0(5)=[];g0(3)=[];g0(2)=[];g0(1)=[];
        alt_clust(9)=[];alt_clust(6)=[];alt_clust(5)=[];
        alt_clust(3)=[];alt_clust(2)=[];alt_clust(1)=[];
        
        %Diagonalize hessian, replace D with |D|
        [V,D]=eig(H);
        H=V*abs(D)*V^(-1);
        try
            step=-H\g0';
        catch
            return
        end
        
        if norm(step)<tol
            return
        end
        
        alt_clust=alt_clust+step';
        
        %Reformat cluster with appropriate zeros in 1,2,3,5,6,9.
        cluster=[0 0 0 alt_clust(1) 0 0 alt_clust(2:3) 0 alt_clust(4:length(alt_clust))];
        counter=counter+1;
    end
end