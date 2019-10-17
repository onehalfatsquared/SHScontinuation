%Scale down rho from 50 to 1 in steps of alpha (0.1) and record the clusters
%every integer value of rho 
clear

%Read data from file
n=9;                                    %Number of particles
filename=strcat('n',num2str(n),'adjust.txt'); %files are nx.txt s.t. x= num part.
fileID=fopen(filename,'r');             %Open file
formatSpec='%f';                        %Format as float
data=fscanf(fileID,formatSpec);         %read in data    
fclose(fileID);                         %Close file

%Seperate column vector into an array of clusters, num clusters by 3n.
num_clusters=length(data)/(3*n);
clusters=reshape(data,[3*n,num_clusters])';

rho=50; E=0.1; rho_end=1;
for alpha=0.1:0.1:rho-rho_end
    for i=1:10%num_clusters
        cluster=clusters(i,:);       %Cluster interested in
        err_catch=cluster;           %Store original cluster in case error
        initE=MP(c2p(cluster),rho-alpha,E); %Energy of start config in new potential
        
        %Perform preconditioned cg step, check if lower energy
        [cluster,flag]=preConjGradDesc(cluster,rho-alpha,E);
        Ediff=MP(c2p(cluster),rho-alpha,E)-initE;
        if Ediff>0 || isnan(Ediff) || flag==1
            %If pCG doesnt work, reset and use GD
            cluster=err_catch;
            cluster=gradDesc(cluster,c2p(cluster),rho-alpha,E);
            Ediff=MP(c2p(cluster),rho-alpha,E)-initE;
            if Ediff>0 || isnan(Ediff)
                %If GD also doesn't work, use starting coords
                cluster=err_catch;
                fprintf('Cluster %d ran into an error at rho=%.2f. Logged starting coords.\n',i,rho-alpha);
            end
        end
        
        %Check for saddle points, get to min, update cluster
        eigvals=eig(hessMorse(cluster,rho-alpha,E));
        negs=length(find(eigvals<0));
        if negs>0
            fprintf('Cluster %d has %d negative eigenvalue(s).\n',i,negs)
            try
                err_catch=cluster;
                cluster=escapeSaddle(cluster,rho-alpha,E);
                cluster=sfNewton(cluster,rho-alpha,E);
            catch
                cluster=err_catch;
                fprintf('Couldnt escape saddle\n')
            end
        end
        clusters(i,:)=cluster;
    end
    if rem(alpha,1)<1e-2
        %rho is an integer, write clusters to a file
        rho-alpha
        filename=strcat('n',num2str(n),'rho',num2str(round(rho-alpha)),'NM.txt'); %File name nxrhoy.txt
        fileID=fopen(filename,'w');                     %Open file
        for i=1:num_clusters
            fprintf(fileID,'%.16f \n',clusters(i,:));      %Write cluster
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Local Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function [cluster, counter]=conjGradDesc(cluster, rho, E)
    %Perform conjugate gradient descent 
    cutoff=250;   %Cut off if doesn't converge
    counter=0;    %Count number of iterations
    tol=1e-15;    %Tolerance for updates
    
    %Set initial search direction as negative of gradient
    g0=-grad(cluster,rho,E); %Gradient of potential
    h=g0;        %First search direction
    
    
    while counter<cutoff
        sigma=tol/max(abs(h)); %adaptive finite diff
        %Perform line min to find optimal step size. Secant method.
        gsig=grad(cluster+sigma*h,rho,E);
        lambda=sigma*g0*h'/(gsig*h'+g0*h');
        if lambda==Inf
            return
        end
        
        %Stop if tol is met. stepsize*dir<tol
        mag=norm(h);
        if mag*lambda<tol
            return
        end
        
        %Update position of the min
        cluster=cluster+lambda*h;
        
        %Compute gradient locally at new point
        g1=-grad(cluster,rho,E);
        
        %Compute gram-schmidt coeff and update search dir
        %gamma=(g1*g1')/(g0*g0');       %Flecther-Reeves update
        gamma=((g1-g0)*g1')/(g0*g0');  %Polak-Ribierre update
        h=g1+gamma*h;                  %Conjugate search direction
        g0=g1;                         %Update gradient at previous
        counter=counter+1;             %Increment step counter
    end
end

function clust=escapeSaddle(clust, rho,E)
    %Escape a saddle point by searching in direction of negative eigenvalue
    %and performing CG step repeatedly.
    counter=0; n=length(clust); max_iters=100;
    while counter<max_iters
        H=hessMorse(clust,rho,E); %Hessian
        [V,D]=eig(H);           %Eigenvalue decomp
        nr=ceil(find(D<0)/(n)); %index of relevant eigenvector
        if isempty(nr)
            break               %End if negative eigenvalues are gone
        end
        ind=ceil(find(D(nr,nr)==max(diag(D(nr,nr))))/length(nr)); %smallest eigvec
        %ind=ceil(find(D(nr,nr)==min(diag(D(nr,nr))))/length(nr));  %Largest eigvec
        %ind=length(nr)
        nr=nr(ind);             %Start with smallest eigenvector in case there are several
        eigvec=V(:,nr)';        %Eigenvector corresponding to neg eigval
       
        
        %line search - newton method
        g0=-grad(clust,rho,E); %Gradient of potential
        delta=-g0*eigvec'/(eigvec*H*eigvec');
        clust=clust+delta*eigvec;
        
        %Grad desc
        %[clust]=conjGradDesc(clust,rho,E);
        clust=gradDesc(clust,c2p(clust),rho,E);
        counter=counter+1;
    end
    [clust]=conjGradDesc(clust,rho,E);
    if counter==max_iters
        fprintf('Could not escape saddle\n')
    end
end

function [cluster,flag]=preConjGradDesc(cluster, rho, E)
    %Perform conjugate gradient descent with preconditioner
    inner=100;   %Inner cut off if doesn't converge
    outer=50;   %Outer cutoff if doesn't converge
    counter_i=0;    %Count number of inner iterations
    counter_o=0;    %Count number of outer iterations
    tol=1e-15;    %Tolerance for updates
    d0=0;      %to determine convergence
    orig=cluster; %original cluster
    flag=0;
    
    while counter_o<outer
        %Compute Hessian, create preconditioner on nonzero eig subspace
        H=hessMorse(cluster,rho,E);
        
        %Delete zero rows and columns 1,2,3,5,6,9. Compute ichol
        H( ~any(H,2), : ) = [];  %rows
        H( :, ~any(H,1) ) = [];  %columns
        H=sparse(H); 
        try
            L=ichol(H,struct('type','ict','droptol',1e-3,'diagcomp',1));
            M=L*L';
        catch
            flag=1;
            return 
        end
            
        
        %Set initial search direction as negative of gradient times precond.
        %Delete the zero entries 1,2,3,5,6,9
        g0=-grad(cluster,rho,E); %Gradient of potential
        alt_clust=cluster;
        g0(9)=[];g0(6)=[];g0(5)=[];g0(3)=[];g0(2)=[];g0(1)=[];
        alt_clust(9)=[];alt_clust(6)=[];alt_clust(5)=[];
        alt_clust(3)=[];alt_clust(2)=[];alt_clust(1)=[];
        g0tilde=M\g0'; g0tilde=g0tilde';
        h0=g0tilde;        %First search direction
        
        while counter_i<inner
            %Compute step
            step=(g0*g0tilde')/(h0*H*h0');
            
            %Stop if tol is met. stepsize*dir<tol
            mag=norm(h0);
            if abs(mag*step)<tol
                break
            end
            
            %Take step
            alt_clust=alt_clust+step*h0;
            
            %Update residual (g)
            g1=g0-step*(H*h0')';
            g1tilde=M\g1'; g1tilde=g1tilde';
            
            %Compute correction factor and update
            beta=(g1*g1tilde')/(g0*g0tilde');
            g0=g1; g0tilde=g1tilde;
            h0=g1tilde+beta*h0;
            counter_i=counter_i+1;
        end
        
        %Reformat cluster with appropriate zeros in 1,2,3,5,6,9.
        cluster=[0 0 0 alt_clust(1) 0 0 alt_clust(2:3) 0 alt_clust(4:length(alt_clust))];
        counter_o=counter_o+1;
        counter_i=0;
        
        %Check for convergence
        d1=MP(c2p(cluster),rho,E)-MP(c2p(orig),rho,E);
        if abs(d1-d0)<tol
            break
        else
            d0=d1;
        end
    end
end

function [cluster, counter]=conjGradDescRes(cluster, rho, E)
    %Perform conjugate gradient descent with restarts
    comps=length(cluster)-6;     %Numer of components
    cutoff=10*comps;   %Cut off if doesn't converge
    counter=0;     %Count number of iterations
    tol=1e-14;    %Tolerance for updates
    
    while counter<cutoff
        %Set initial search direction as negative of gradient
        g0=-grad(cluster,rho,E); %Gradient of potential
        h=g0;        %First search direction
        resetcounter=0; %Count number of conjugate directions used
        while resetcounter<comps
            sigma=tol/max(abs(h)); %adaptive finite diff
            %Perform line min to find optimal step size. Secant method.
            gsig=grad(cluster+sigma*h,rho,E);
            
            lambda=sigma*g0*h'/(gsig*h'+g0*h');
            if lambda==Inf
                return
            end
            
            %Stop if tol is met. stepsize*dir<tol
            mag=norm(h);
            if mag*lambda<tol
                return
            end
            
            %Update position of the min
            cluster=cluster+lambda*h;
            
            %Compute gradient locally at new point
            g1=-grad(cluster,rho,E);
            
            %Compute gram-schmidt coeff and update search dir
            %gamma=(g1*g1')/(g0*g0');       %Flecther-Reeves update
            gamma=((g1-g0)*g1')/(g0*g0');  %Polak-Ribierre update
            h=g1+gamma*h;                  %Conjugate search direction
            g0=g1;                         %Update gradient at previous
            counter=counter+1;             %Increment step counter
            resetcounter=resetcounter+1;   %Update reset counter
        end
    end
end



function [cluster,particles,counter]=gradDesc(cluster, particles, rho, E)
    %Perform gradient descent 
    cutoff=750;                      %Cut off if doesn't converge
    counter=0;                       %Count number of iterations
    tol=1e-15;                       %tol to end 
    %perform the first gradient descent step
    g0=grad(cluster,rho,E); %Gradient of potential
    alpha=0.85;                       %Parameter to exp. for step
        for i=0:100
            L=MP(c2p(cluster-alpha^i*g0),rho,E);
            R=MP(particles,rho,E);
            if L<=R
                cluster1=cluster-alpha^i*g0; %first update to cluster
                particles=c2p(cluster1);    %first update to particles
                break
            elseif i==100 
                %Could not find step size to make U smaller
                return
            end
        end
    counter=1;
    while counter<cutoff
        g1=grad(cluster1,rho,E); %Gradient of potential 
        %determine step size and update 
        dx=cluster1-cluster;
        dg=g1-g0;
        lambda=(dg*dx')/(dg*dg');         %Step size
        mag=norm(g1);                      %magnitude of gradient
        cluster=cluster1; g0=g1;           %Update k-1 elements
        cluster1=cluster1-lambda*g1;       %Update k elements
        particles=c2p(cluster1);          %Update particles list
        counter=counter+1;                 %Increment counter
        if mag*lambda<tol
            return
        end
    end
end