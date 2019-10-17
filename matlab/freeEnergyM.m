function F = freeEnergyM(clusts,clusters,rho,E,T)
    %Evaluate the free energy of a Morse cluster. 
    %clusts contains index of clusters to evaluate FE for
    %clusters is array of all clusters for given N.
    %rho, E are parameters to morse potential. T is temp. 
    
    %parameters
    kb=1.38064852*10-23; %boltzmann constant
    beta = 1;
    
    %first compute all probs and normalizing factor
    [num_clust,~]=size(clusters); %number of clusters
    p=zeros(num_clust,1); %initialize probability storage
    for i=1:num_clust
        p(i)=computeProb(clusters(i,:),rho,E,T);
    end
    Z = sum(p);
    
    %compute the free energies and return values
    f=-1/beta*log(p); 
    F=f(clusts);  
end

function P = computeProb(clust,rho,E,T)
    %Compute probability up to constants of a cluster
    
    kb=1.38064852*10-23; %boltzmann constant
    beta=1/(kb*T); %inverse temp
    beta = 1;
    Zr = rotPF(clust);
    Zv = vibPF(clust,rho,E);
    Z=Zr*Zv; %partition function up to constants
    U=MP(c2p(clust),100,E); %morse potential energy
    boltz=exp(-beta*U); %boltzmann factor
    P=Z*boltz; %probability, un-normalized
end

function zR = rotPF(clust)
    %evaluate the rotational partition function of a cluster
    
    M=inertiaTensor(clust); I=det(M); %inertia tensor + determinant
    s=symNum(clust); %evaluate symmetry number
    zR=factorial(7)*sqrt(I)/s;  %eval partition fn
end

function zV = vibPF(clust,rho,E)
    %evaluate vibrational partition function
    
    n=length(clust)-6; %num of Dofs
    H=hessMorse(clust,rho,E);  %hessian with zeros eigs
    e=eig(H); e=e(e>0); e=prod(e); %remove zero eigs, take product
    zV=(2*pi)^(n/2)/sqrt(e); %evaluate vib partition fn
end

function s = symNum(clust)
    %evaluate the symmetry number of a cluster
    
    tol=1e-5; %tolerance for same matrix elements
    particles=c2p(clust); parts=length(particles); % particle array and num particles
    ipdMatrix = pDists(particles); % inter particle distance matrix
    AP=perms(1:parts); %list of all perms
    s=0; %initialize sym number
    %Loop to check how many perms give same ipd matrix
    for i=1:factorial(parts)
        test=particles(AP(i,:),:);
        testMatrix=pDists(test); 
        if max(max(abs(testMatrix-ipdMatrix)))<tol
            s=s+1;
        end
    end
end