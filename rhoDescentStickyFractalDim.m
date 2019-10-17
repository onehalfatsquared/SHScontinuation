%Scale down rho from 50 to 1 in steps of alpha (0.1). Estimate fractal dimension of path
%Keep constant sticky parameter. 
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

%Define parameters to descent, initial energy, and initial sticky param.
rho=50; rho_end=1; eigTol=-1e-10; Etol=1e-14; %Define parameters
E=7.5; %Set the initial energy. Must be >5.7.  
kap=stickyEval(E,rho); %Evaluate the initial sticky parameter
i=1; %Cluster interested in 
data=[]; %initialize data storage. record cluster at every rho value
counter=1;

%Perform rho descent. 
for alpha=0:0.1:rho-rho_end
    %Determine the energy  for this rho value to keep kappa const.
    E=stickyNewton(E, rho-alpha, kap);
    
    cluster=clusters(i,:);       %Cluster interested in
    err_catch=cluster;           %Store original cluster in case error
    
    %Perform cg step, check if lower energy
    [cluster,c]=conjGradDesc(cluster,rho-alpha,E);
    if errorCheck(cluster,err_catch,rho-alpha,E)==0
        %If CG doesnt work, reset cluster and call the errorDescent
        cluster=err_catch;
        [cluster,f]=errorDescent(cluster,rho-alpha,E);
        if f==1
            %If errorDescent doesnt work, use starting coords
            fprintf('Cluster %d encountered an error at rho=%.2f. Logged starting coords.\n',i,rho-alpha);
        end
    end
    
    %Check for saddle points, get to min, update cluster
    eigvals=eig(hessMorse(cluster,rho-alpha,E));
    negs=length(find(eigvals<eigTol));
    if negs>0
        fprintf('Cluster %d has %d negative eigenvalue(s) at rho=%.2f.\n',i,negs,rho-alpha)
        norm(grad(cluster,rho-alpha,E))
        cluster=escapeSaddle(cluster,rho-alpha,E);
        if norm(cluster)>10
            cluster=err_catch;
        end
    end
    clusters(i,:)=cluster;
    data(counter,:)=cluster; counter=counter+1;
    %norm(grad(cluster,rho-alpha,E))
end
fprintf('Optimization concluded. Computing correlation dimension.\n')

%points are stored in data. estimate fractal dimension by correlation
%dimension. 
N=size(data,1); %Number of points along path
r=logspace(-10,0,11); %log spaced vector 
C=zeros(1,length(r)); %initialize corr fn
%Loop over all data pairwise
for a=1:N
    for b=1:N
        if a~=b
            d=norm(data(a,:)-data(b,:));
            if d>1e-14
                theta=(r-d)>0;
                C=C+theta; 
            end
        end
    end
end
C=C/(N*(N-1));
loglog(r,C)
p=polyfit(log(r),log(C),1)
            



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Local Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cluster, counter]=conjGradDesc(cluster, rho, E)
    %Perform conjugate gradient descent 
    cutoff=2100;   %Cut off if doesn't converge
    counter=0;    %Count number of iterations
    tol=1e-14;    %Tolerance for updates on gradient
    
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
        mag=norm(g0);
        if mag<tol
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

function [clust_out, counter]=conjGradDescRes(cluster, rho, E)
    %Perform conjugate gradient descent with restarts
    comps=length(cluster)-6;     %Numer of components
    cutoff=1000;   %Cut off if doesn't converge
    counter=0;     %Count number of iterations
    gtol=1e-14;    %gradient tolerance for updates
    eTol=0;   %energy tolerance for updates
    e1=MP(c2p(cluster),rho,E);
    
    while counter<cutoff
        %Set initial search direction as negative of gradient
        g0=-grad(cluster,rho,E); %Gradient of potential
        h=g0;        %First search direction
        resetcounter=0; %Count number of conjugate directions used
        clust_out=cluster; %set cluster to output
        while resetcounter<comps
            sigma=gtol/max(abs(h)); %adaptive finite diff
            %Perform line min to find optimal step size. Secant method.
            gsig=grad(cluster+sigma*h,rho,E);
            
            lambda=sigma*g0*h'/(gsig*h'+g0*h');
            if lambda==Inf
                return
            end
            
            %Stop if tol is met. stepsize*dir<tol
            mag=norm(h);
            if mag<gtol
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
        e2=MP(c2p(cluster),rho,E); ediff=e2-e1; e1=e2;
        if ediff>eTol
            return
        end
    end
end

function [cluster,particles,counter]=gradDesc(cluster, particles, rho, E)
    %Perform gradient descent 
    cutoff=1000;                      %Cut off if doesn't converge
    counter=0;                       %Count number of iterations
    tol=1e-10;                       %tol to end 
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
        if mag<tol
            return
        end
    end
end

function left=escapeSaddle(clust, rho,E)
    %Escape a saddle point by searching in direction of negative eigenvalue
    %and performing CG step. 
    eigTol=-1e-7; initialStep=1e-10; 
    
    H=hessMorse(clust,rho,E); %Hessian matrix
    [V,D]=eig(H); D=diag(D);  %Eigenvalue decomposition
    nr=ceil(find(D<eigTol));  %index of relevant eigenvectors
    D(nr)
    if length(nr)==1
        %One negative eigenvalue. Search both directions. 
        eigvec=V(:,nr)';      %Eigenvector corresponding to neg eigval
        
        %Search left
        step=5*initialStep;
        while step<1
            left=clust-step*eigvec;  %step along eigvec
            left=conjGradDesc(left,rho,E); if isnan(left(1)); left=clust; end  %perform minimization
            Hleft=hessMorse(left,rho,E); eleft=eig(Hleft);
            if min(eleft)<eigTol || norm(left)>7
                step=step*10
            else
                break
            end
        end
        if step>=1
            fprintf('Could not escape left')
            left=clust;
        end
        
        %Search right
        step=5*initialStep;
        while step<1
            right=clust+step*eigvec;  %step along eigvec
            right=conjGradDesc(right,rho,E);  if isnan(right(1)); right=clust; end  %perform minimization
            Hright=hessMorse(right,rho,E); eright=eig(Hright);
            if min(eright)<eigTol || norm(right)>7
                step=step*10
            else
                break
            end
        end
        if step>=1
            fprintf('Could not escape right')
        end
    end
    
    %Test if the two clusters are the same
    left, right
    testSame(left,right)
    
end

function [clust,flag]=errorDescent(clust, rho,E)
    %If CG runs into an error, perform various attempts at reaching a minimum
    %Try initial grad desc, permuting particles, perturbations, 
    
    err_catch=clust; %save the original cluster in case of error
    n=length(clust)/3; %Length of cluster
    flag=0;
    
    
    %Method 1: GD then CG
    clust=gradDesc(clust,c2p(clust),rho,E); clust=conjGradDesc(clust,rho,E); 
    if errorCheck(clust,err_catch,rho,E)==1
        fprintf('Success at 1\n')
        return
    else
        clust=err_catch;
    end
    
    %Method 2: CG with restarts then CG
    clust=conjGradDescRes(clust,rho,E); clust=conjGradDesc(clust,rho,E); 
    if errorCheck(clust,err_catch,rho,E)==1
        fprintf('Success at 2\n')
        return
    else
        clust=err_catch;
    end
    
    %Method 3: Permute particles, randomly among the last n-3. Try 5 swaps.
    for attempt=1:5
        %Pick index of particles to swap, 0=n, 1=n-1, etc...
        r=randperm(n-3,n-3)-1; x=r(1); y=r(2);  
        %swap particle n-x with particle n-y
        temp=clust(3*n-(2+3*x):3*n-3*x); 
        clust(3*n-(2+3*x):3*n-3*x)=clust(3*n-(2+3*y):3*n-3*y); 
        clust(3*n-(2+3*y):3*n-3*y)=temp; 
        clust=conjGradDesc(clust,rho,E); 
        if errorCheck(clust,err_catch,rho,E)==1
            temp=clust(3*n-(2+3*x):3*n-3*x);
            clust(3*n-(2+3*x):3*n-3*x)=clust(3*n-(2+3*y):3*n-3*y);
            clust(3*n-(2+3*y):3*n-3*y)=temp;
            fprintf('Success at 3\n')
            return
        else
            clust=err_catch;
        end
    end
    
    %Method 4: 10 Random perturbation to all coordinates but 1,2,3,5,6,9. 
    for attempt=1:10
        p=[0,0,0,rand()*1e-12,0,0,rand()*1e-12, rand()*1e-12,0];
        p=[p rand(1,3*n-9)*1e-12]; clust=clust+p;
        clust=conjGradDesc(clust,rho,E); 
        if errorCheck(clust,err_catch,rho,E)==1
            fprintf('Success at 4\n')
            return
        else
            clust=err_catch;
        end
    end
    
    %If all methods fail, return starting cluster, set fail flag to 1
    clust=err_catch; flag=1;
    return 
end


function ok=errorCheck(clust, orig, rho,E)
    %Check if optimization was succesful
    Etol=1e-14; normTol=7;   %Parameters to see if energy decreased, and tolerance for CG blowup
    ok=1;  %initialize output as error free. Set to 0 if error. 
    
    Einit=MP(c2p(orig),rho,E); Efinal= MP(c2p(clust),rho,E); %Initial and final energy
    Ediff=Efinal-Einit;   %energy difference
    
    if Ediff>Etol || isnan(Ediff) || norm(clust)>normTol
        %Ediff, norm(clust), norm(grad(clust,rho,E))
        ok=0;
    end
end




