%get data for jump scatter plot - with rmin and cluster index
clear

%set parameters
n = 9; %done up to n=9, do not redo

%load the file with rho marge values and jump
filename = strcat('n',num2str(n),'fM2.mat'); %File name nxrhoy.txt
load(filename)
load('jumpSP6')


%get number of clusters
[L,~] = size(mg)

%initialize Jump vector
% jumpSP5{1,1} = [];
% for i=1:50
%     jumpSP5{i,1} = [];
% end


% %get the data
% for i =1:L
%     i
%     %get the position of the cluster at rho value just before merge and
%     %after
%     rho = mg{i,1}+1;
%     jumps=[];
%     for clusts = mg{i,2}
%         cb = getCluster(n,rho,clusts); 
%         ca = getCluster(n,rho-1,clusts);
%         [~,~,~,~,~,d] = testSame(ca,cb,1);
%         jumps = [jumps d];
%     end
%     jumpSP5{end+1,1} = [jumps];
%     jumpSP5{end,2} = rho -1;
%     
%    
% end


%get the data
for i =1:L
    i
    %get the position of the cluster at rho value just before merge
    rho = mg{i,1}+1;
    if rho >35
        jumps=[];
        for clusts = mg{i,2}
            cb = getCluster(n,rho,clusts); 
            jumpList = rhoDesc(cb, rho);
            d = max(jumpList);
            jumps = [jumps d];
        end
        jumpSP6{end+1,1} = [jumps];
        jumpSP6{end,2} = rho -1; 
    end
end

%save('jumpSP5.mat','jumpSP5')

function [jumpList] = rhoDesc(cb, rho)
%perform descent
%Define parameters to descent, initial energy, and initial sticky param.
rho_end=rho-1; eigTol=-1e-10; 
energies=[7.5,8.3,9.05]; 
test=1; %1, 2, or 3 correspond to low, med, high energy test
E=energies(test); %Set the initial energy. Must be >5.7.
kap=stickyEval(E,rho); %Evaluate the initial sticky parameter
fails=0;
before = cb; 
jumpList=[];

if rho>10 
    step = 0.1;
else
    step = 0.01;
end

%Perform rho descent.
for alpha=step:step:rho-rho_end
    %Determine the energy  for this rho value to keep kappa const.
    E=stickyNewton(E, rho-alpha, kap);
    err_catch=before;           %Store original cluster in case error
    
    %Perform cg step, check if lower energy
    [after]=conjGradDesc(cb,rho-alpha,E);
    if errorCheck(after,err_catch,rho-alpha,E)==0
        %If CG doesnt work, reset cluster and call the errorDescent
        after=err_catch;
        [after,f]=errorDescent(after,rho-alpha,E);
        if f==1
            %If errorDescent doesnt work, use starting coords
            fprintf('Cluster %d encountered an error at rho=%.2f. Logged starting coords.\n',0,rho-alpha);
            fails=fails+1;
        end
    end
    
    %Check for saddle points, get to min, update cluster
    eigvals=eig(hessMorse(after,rho-alpha,E));
    negs=length(find(eigvals<eigTol));
    if negs>0
        fprintf('Cluster %d has %d negative eigenvalue(s) at rho=%.2f.\n',0,negs,rho-alpha)
        norm(grad(after,rho-alpha,E))
        after=escapeSaddle(after,rho-alpha,E);
        if norm(after)>10
            after=err_catch;
        end
    end
    
    %get the rmsd bewtween before and after
    [~,~,~,~,~,d] = testSame(before,after,1);
    jumpList = [jumpList d];
    before = after; 
end
end


function rc = rmin(cluster, rho)
    %get rmin value
    comp_val = [1.01, 1.06];
    if rho>30
        comp = 1;
    else
        comp = 2;
    end
    
    di=unique(pDists(c2p(cluster))); di=di(di>comp_val(comp)); %use 1.01 for rho 38, 1.05 for rho 6
    rc = min(di);
end

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

function right=escapeSaddle(clust, rho,E)
    %Escape a saddle point by searching in direction of negative eigenvalue
    %and performing CG step. 
    eigTol=-1e-10; initialStep=1e-10; 
    
    H=hessMorse(clust,rho,E); %Hessian matrix
    [V,D]=eig(H); D=diag(D);  %Eigenvalue decomposition
    nr=ceil(find(D<eigTol));  %index of relevant eigenvectors
    D(nr)
    if length(nr)==1
        %One negative eigenvalue. Search both directions. 
        eigvec=V(:,nr)';      %Eigenvector corresponding to neg eigval
    elseif length(nr)>1
        ind=find(D==min(D));  %Max eigvec
        %ind=find(round(D,12)==round(1/min(1./D),12)); %Min eigvec
        eigvec=V(:,ind)';  %#ok<FNDSB>
    end
        
    %Search left
    step=5*initialStep;
    while step<1
        left=clust-step*eigvec;  %step along eigvec
        left=errorDescent(left,rho,E); if isnan(left(1)); left=clust; end  %perform minimization
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
        right=errorDescent(right,rho,E);  if isnan(right(1)); right=clust; end  %perform minimization
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
    Etol=1e-14; normTol=6;   %Parameters to see if energy decreased, and tolerance for CG blowup
    ok=1;  %initialize output as error free. Set to 0 if error. 
    
    Einit=MP(c2p(orig),rho,E); Efinal= MP(c2p(clust),rho,E); %Initial and final energy
    Ediff=Efinal-Einit;   %energy difference
    
    if Ediff>Etol || isnan(Ediff) || norm(clust)>normTol
        %Ediff, norm(clust), norm(grad(clust,rho,E))
        ok=0;
    end
end

